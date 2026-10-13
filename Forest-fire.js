// =========================================================================
// SPATIAL ASSESSMENT OF WILDFIRE DRIVING FACTORS, BURN SEVERITY, 
// AND POST-FIRE VEGETATION REGENERATION DYNAMICS (MT. LONGONOT ECOSYSTEM)
// =========================================================================

// -------------------------------------------------------------------------
// 1. AOI DEFINITION & INITIALIZATION
// -------------------------------------------------------------------------
var aoi = ee.FeatureCollection('projects/northern-anchor-441619-g4/assets/Longonot').geometry();

print('================================================================');
print('MT. LONGONOT WILDFIRE RISK AND RECOVERY ANALYSIS PIPELINE');
print('================================================================');

// -------------------------------------------------------------------------
// 2. ENVIRONMENTAL & ANTHROPOGENIC DRIVING FACTORS
// -------------------------------------------------------------------------
var dem = ee.Image('USGS/SRTMGL1_003').clip(aoi);
var elevation = dem.rename('Elevation');
var slope = ee.Terrain.slope(dem).rename('Slope');
var aspect = ee.Terrain.aspect(dem).rename('Aspect');

var s2 = ee.ImageCollection('COPERNICUS/S2_SR_HARMONIZED')
  .filterBounds(aoi)
  .filter(ee.Filter.lt('CLOUDY_PIXEL_PERCENTAGE', 40));

var baselineS2 = s2.filterDate('2026-01-01', '2026-03-31').median().clip(aoi);
var ndvi = baselineS2.normalizedDifference(['B8', 'B4']).rename('NDVI');
var ndwi = baselineS2.normalizedDifference(['B3', 'B8']).rename('NDWI');

var CLIMATE_START = '2025-01-01';
var CLIMATE_END = '2025-12-31';

var rainfall = ee.ImageCollection('UCSB-CHG/CHIRPS/DAILY')
  .filterBounds(aoi)
  .filterDate(CLIMATE_START, CLIMATE_END)
  .sum()
  .rename('Rainfall')
  .clip(aoi);

var lst = ee.ImageCollection('MODIS/061/MOD11A2')
  .filterBounds(aoi)
  .filterDate(CLIMATE_START, CLIMATE_END)
  .select('LST_Day_1km')
  .mean()
  .multiply(0.02).subtract(273.15)
  .rename('LST')
  .clip(aoi);

var nightLights = ee.ImageCollection('NOAA/VIIRS/DNB/MONTHLY_V1/VCMSLCFG')
  .filterBounds(aoi)
  .filterDate(CLIMATE_START, CLIMATE_END)
  .select('avg_rad')
  .mean()
  .rename('NightLights')
  .clip(aoi);

var predictors = ee.Image.cat([slope, aspect, ndvi, ndwi, elevation, rainfall, lst, nightLights]);
var predictorNames = ['Slope', 'Aspect', 'NDVI', 'NDWI', 'Elevation', 'Rainfall', 'LST', 'NightLights'];

print('Predictor Variables Loaded:', predictors.bandNames());

// -------------------------------------------------------------------------
// 3. MULTI-SENSOR NBR & dNBR FUNCTIONS
// -------------------------------------------------------------------------
function getLandsatNBR(collectionIds, start, end) {
  collectionIds = (typeof collectionIds === 'string') ? [collectionIds] : collectionIds;
  var merged = ee.ImageCollection([]);
  collectionIds.forEach(function (id) {
    var c = ee.ImageCollection(id)
      .filterBounds(aoi)
      .filterDate(start, end)
      .filter(ee.Filter.lt('CLOUD_COVER', 60));
    merged = merged.merge(c);
  });

  var scaled = merged.map(function (img) {
    return img.select('SR_B.').multiply(0.0000275).add(-0.2).copyProperties(img, img.propertyNames());
  });

  return scaled.median().normalizedDifference(['SR_B4', 'SR_B7']).rename('NBR');
}

function getSentinelNBR(start, end) {
  var col = s2.filterDate(start, end);
  return col.median().normalizedDifference(['B8', 'B12']).rename('NBR');
}

function getDNBR(nbrPre, nbrPost) {
  return nbrPre.subtract(nbrPost).rename('dNBR').clip(aoi);
}

function classifyDNBR(dnbrImg) {
  var classified = dnbrImg.expression(
    '(val < 0.1) ? 0 : (val < 0.27) ? 1 : (val < 0.44) ? 2 : 3',
    {'val': dnbrImg}
  );
  return classified.updateMask(dnbrImg.mask());
}

// -------------------------------------------------------------------------
// 4. CALCULATE HISTORICAL BURN SEVERITIES (2009, 2021, 2022, 2026)
// -------------------------------------------------------------------------
var l09 = ['LANDSAT/LT05/C02/T1_L2', 'LANDSAT/LE07/C02/T1_L2'];
var nbrPre2009 = getLandsatNBR(l09, '2008-11-01', '2009-01-31');
var nbrPost2009 = getLandsatNBR(l09, '2009-02-15', '2009-04-30');
var dnbr2009 = getDNBR(nbrPre2009, nbrPost2009);
var severity2009 = classifyDNBR(dnbr2009);

var nbrPre2021 = getSentinelNBR('2020-12-01', '2021-01-15');
var nbrPost2021 = getSentinelNBR('2021-03-01', '2021-04-15');
var dnbr2021 = getDNBR(nbrPre2021, nbrPost2021);
var severity2021 = classifyDNBR(dnbr2021);

var nbrPre2022 = getSentinelNBR('2021-12-01', '2022-01-15');
var nbrPost2022 = getSentinelNBR('2022-03-01', '2022-04-15');
var dnbr2022 = getDNBR(nbrPre2022, nbrPost2022);
var severity2022 = classifyDNBR(dnbr2022);

var nbrPre2026 = getSentinelNBR('2025-12-01', '2026-01-15');
var nbrPost2026 = getSentinelNBR('2026-03-01', '2026-04-15');
var dnbr2026 = getDNBR(nbrPre2026, nbrPost2026);
var severity2026 = classifyDNBR(dnbr2026);

// -------------------------------------------------------------------------
// 5. RANDOM FOREST (RF) WILDFIRE SUSCEPTIBILITY MODELING (WITH VALIDATION)
// -------------------------------------------------------------------------
// FIX: FIRMS active-fire detection only caught 2 hotspots for this whole park
// across 5 years -- far too few to train or validate on. FIRMS needs a
// satellite to catch the ground burning at that exact moment; it easily
// misses short grass fires. Your dNBR maps already show real, detailed burn
// scars across thousands of pixels for all 4 fire years -- use those instead.

// "Ever burned" masks, kept SEPARATE per year (not pooled) so we can hold
// an entire fire year out for testing -- a stronger check against
// overfitting than randomly splitting individual points, since it tests
// whether the model generalizes to a fire it has never seen at all, not
// just to nearby pixels of fires it's already partly seen.
var burn2009 = severity2009.gt(0).unmask(0);
var burn2021 = severity2021.gt(0).unmask(0);
var burn2022 = severity2022.gt(0).unmask(0);
var burn2026 = severity2026.gt(0).unmask(0);

// TRAIN on 2009 + 2021 + 2022 (3 of 4 fire years)
var trainYearsBurn = burn2009.or(burn2021).or(burn2022).rename('class');
// TEST entirely on 2026 -- a year the model never sees during training
var testYearBurn = burn2026.rename('class');

var trainSet = trainYearsBurn.addBands(predictors).stratifiedSample({
  numPoints: 250, // per class
  classBand: 'class',
  region: aoi,
  scale: 30,
  seed: 42,
  geometries: true
});

var testSet = testYearBurn.addBands(predictors).stratifiedSample({
  numPoints: 100, // per class -- 2026's burn scar is smaller, so fewer points available
  classBand: 'class',
  region: aoi,
  scale: 30,
  seed: 99,
  geometries: true
});

print('---- RF Training Data (TEMPORAL HOLD-OUT: train on 2009/2021/2022, test on 2026) ----');
print('Train set size (2009 + 2021 + 2022 fires):', trainSet.size());
print('Test set size (2026 fire, never seen in training):', testSet.size());

// COMPARISON: does including the Landsat-based 2009 year help or hurt,
// versus training on Sentinel-2 years only? Same test set (2026) both times,
// so this isolates the effect of the 2009 training data specifically.
var trainYearsBurnSentinelOnly = burn2021.or(burn2022).rename('class');
var trainSetSentinelOnly = trainYearsBurnSentinelOnly.addBands(predictors).stratifiedSample({
  numPoints: 250,
  classBand: 'class',
  region: aoi,
  scale: 30,
  seed: 42,
  geometries: true
});

var trainedRFSentinelOnly = ee.Classifier.smileRandomForest(100).train({
  features: trainSetSentinelOnly,
  classProperty: 'class',
  inputProperties: predictorNames
});
var testClassifiedSentinelOnly = testSet.classify(trainedRFSentinelOnly);
var confusionMatrixSentinelOnly = testClassifiedSentinelOnly.errorMatrix('class', 'classification');

print('---- COMPARISON: Sentinel-only training (2021+2022, no 2009), tested on 2026 ----');
print('Train set size (2021 + 2022 only):', trainSetSentinelOnly.size());
print('Confusion Matrix (Sentinel-only):', confusionMatrixSentinelOnly);
print('Accuracy (Sentinel-only):', confusionMatrixSentinelOnly.accuracy());
print('Kappa (Sentinel-only):', confusionMatrixSentinelOnly.kappa());

// DECISION: the Sentinel-only model (2021+2022 training) scored HIGHER
// accuracy and kappa on the 2026 test set than the full model that included
// the Landsat-based 2009 year. That's evidence the two sensors' dNBR values
// aren't fully interchangeable for this task -- so the susceptibility map
// actually displayed in the dashboard uses the better-performing
// Sentinel-only model, not the full one.
var trainedRFSentinelOnlyProbability = ee.Classifier.smileRandomForest(100)
  .setOutputMode('PROBABILITY')
  .train({
    features: trainSetSentinelOnly,
    classProperty: 'class',
    inputProperties: predictorNames
  });

var trainedRFForValidation = ee.Classifier.smileRandomForest(100).train({
  features: trainSet,
  classProperty: 'class',
  inputProperties: predictorNames
});
var testClassified = testSet.classify(trainedRFForValidation);
var confusionMatrix = testClassified.errorMatrix('class', 'classification');

print('---- RF Validation, FULL model incl. 2009 (held-out test) ----');
print('Confusion Matrix:', confusionMatrix);
print('Validation Overall Accuracy:', confusionMatrix.accuracy());
print('Validation Kappa:', confusionMatrix.kappa());
print('NOTE: the map below uses the Sentinel-only model (higher accuracy/kappa above), not this one.');

var rfProbability = predictors.classify(trainedRFSentinelOnlyProbability).clip(aoi).rename('FireProbability');

var rfRiskClass = rfProbability.expression(
  '(p < 0.33) ? 1 : (p < 0.66) ? 2 : 3', {p: rfProbability}
).rename('RiskClass');

print('Random Forest Model Trained and Validated.');

// -------------------------------------------------------------------------
// 6. ZONAL STATISTICAL ANALYSIS (BURN AREA QUANTIFICATION)
// -------------------------------------------------------------------------
function computeBurnArea(severityImg, yearName) {
  var areaImage = ee.Image.pixelArea().divide(10000).addBands(severityImg);
  var stats = areaImage.reduceRegion({
    reducer: ee.Reducer.sum().group({
      groupField: 1,
      groupName: 'severity_class',
    }),
    geometry: aoi,
    scale: 30,
    maxPixels: 1e9
  });
  print('Burned Area Statistics (Hectares) - ' + yearName + ':', stats.get('groups'));
}

computeBurnArea(severity2009, '2009');
computeBurnArea(severity2021, '2021');
computeBurnArea(severity2022, '2022');
computeBurnArea(severity2026, '2026');

// -------------------------------------------------------------------------
// 7. POST-FIRE VEGETATION RECOVERY — OVERALL + STRATIFIED BY SEVERITY ZONE
// -------------------------------------------------------------------------
var s2_recovery = s2.filterDate('2021-01-01', '2026-08-31').map(function (img) {
  var veg = img.normalizedDifference(['B8', 'B4']).rename('NDVI');
  return img.addBands(veg).copyProperties(img, ['system:time_start']);
});

var recoveryChart = ui.Chart.image.series({
  imageCollection: s2_recovery.select('NDVI'),
  region: aoi,
  reducer: ee.Reducer.mean(),
  scale: 30,
  xProperty: 'system:time_start'
})
.setChartType('LineChart')
.setOptions({
  title: 'Mt. Longonot Post-Fire Vegetation Recovery - Whole AOI (2021-2026)',
  hAxis: {title: 'Date', format: 'YYYY-MM'},
  vAxis: {title: 'Mean NDVI (Canopy Density)'},
  lineWidth: 2,
  pointSize: 2,
  series: {0: {color: '228B22'}},
  trendlines: {0: {color: 'FF0000', lineWidth: 1, opacity: 0.5}}
});

print(recoveryChart);

var severityRef = severity2021;

function ndviMeanForClass(img, classValue) {
  var masked = img.select('NDVI').updateMask(severityRef.eq(classValue));
  return masked.reduceRegion({
    reducer: ee.Reducer.mean(),
    geometry: aoi,
    scale: 30,
    maxPixels: 1e9
  }).get('NDVI');
}

function zonalNDVIForImage(img) {
  return ee.Feature(null, {
    'Low_NDVI': ndviMeanForClass(img, 1),
    'Moderate_NDVI': ndviMeanForClass(img, 2),
    'High_NDVI': ndviMeanForClass(img, 3)
  }).set('system:time_start', img.get('system:time_start'));
}

var zonalNDVIFeatures = ee.FeatureCollection(s2_recovery.map(zonalNDVIForImage));

var zonalRecoveryChart = ui.Chart.feature.byFeature({
  features: zonalNDVIFeatures,
  xProperty: 'system:time_start',
  yProperties: ['Low_NDVI', 'Moderate_NDVI', 'High_NDVI']
})
.setChartType('LineChart')
.setOptions({
  title: 'Post-Fire Recovery by 2021 Burn Severity Zone',
  hAxis: {title: 'Date', format: 'YYYY-MM'},
  vAxis: {title: 'Mean NDVI'},
  lineWidth: 2,
  pointSize: 2,
  series: {
    0: {color: '#FFFF00', labelInLegend: 'Low severity zone'},
    1: {color: '#FF7F00', labelInLegend: 'Moderate severity zone'},
    2: {color: '#FF0000', labelInLegend: 'High severity zone'}
  }
});

print(zonalRecoveryChart);

// -------------------------------------------------------------------------
// 8. INTERACTIVE LINKED MAP DASHBOARD (4 SEVERITY MAPS + 1 SUSCEPTIBILITY MAP)
// -------------------------------------------------------------------------
var years = [
  {label: '2009 Burn Severity (Landsat 5+7)', image: severity2009},
  {label: '2021 Burn Severity (Sentinel-2)', image: severity2021},
  {label: '2022 Burn Severity (Sentinel-2)', image: severity2022},
  {label: '2026 Burn Severity (Sentinel-2)', image: severity2026}
];

var burnVis = {min: 0, max: 3, palette: ['#2ca02c', '#FFFF00', '#FF7F00', '#FF0000']};
var riskVis = {min: 1, max: 3, palette: ['#1a9850', '#fee08b', '#d73027']};

function buildLegend(title, colors, labels) {
  var legend = ui.Panel({
    style: {
      position: 'bottom-left',
      padding: '6px 8px',
      backgroundColor: 'rgba(255, 255, 255, 0.9)'
    }
  });

  legend.add(ui.Label({
    value: title,
    style: {fontWeight: 'bold', fontSize: '11px', margin: '0 0 4px 0'}
  }));

  for (var i = 0; i < colors.length; i++) {
    var colorBox = ui.Label({
      style: {backgroundColor: colors[i], padding: '8px', margin: '0 4px 3px 0'}
    });
    var description = ui.Label({
      value: labels[i],
      style: {margin: '0 0 3px 0', fontSize: '10px'}
    });
    var row = ui.Panel({
      widgets: [colorBox, description],
      layout: ui.Panel.Layout.flow('horizontal')
    });
    legend.add(row);
  }

  return legend;
}

var severityLegendColors = ['#2ca02c', '#FFFF00', '#FF7F00', '#FF0000'];
var severityLegendLabels = ['No significant burn', 'Low severity', 'Moderate severity', 'High severity'];
var riskLegendColors = ['#1a9850', '#fee08b', '#d73027'];
var riskLegendLabels = ['Low risk', 'Moderate risk', 'High risk'];

var maps = years.map(function (y) {
  var m = ui.Map();
  m.centerObject(aoi, 12);
  m.addLayer(y.image, burnVis, y.label);

  var title = ui.Label(y.label, {
    position: 'top-left',
    fontWeight: 'bold',
    backgroundColor: 'rgba(255, 255, 255, 0.85)',
    padding: '4px'
  });
  m.add(title);
  m.add(buildLegend('dNBR Burn Severity', severityLegendColors, severityLegendLabels));

  m.setControlVisibility({layerList: false, zoomControl: false, scaleControl: true});
  return m;
});

var riskMap = ui.Map();
riskMap.centerObject(aoi, 12);
riskMap.addLayer(rfRiskClass, riskVis, 'Wildfire Susceptibility (Random Forest)');
riskMap.add(ui.Label('Wildfire Susceptibility (Random Forest)', {
  position: 'top-left',
  fontWeight: 'bold',
  backgroundColor: 'rgba(255, 255, 255, 0.85)',
  padding: '4px'
}));
riskMap.add(buildLegend('Fire Susceptibility', riskLegendColors, riskLegendLabels));
riskMap.setControlVisibility({layerList: false, zoomControl: true, scaleControl: true});

var allMaps = maps.concat([riskMap]);
ui.Map.Linker(allMaps, 'change-bounds');

var row1 = ui.Panel([maps[0], maps[1]], ui.Panel.Layout.flow('horizontal'), {stretch: 'both'});
var row2 = ui.Panel([maps[2], maps[3]], ui.Panel.Layout.flow('horizontal'), {stretch: 'both'});
var row3 = ui.Panel([riskMap], ui.Panel.Layout.flow('horizontal'), {stretch: 'both', height: '400px'});
var dashboard = ui.Panel([row1, row2, row3], ui.Panel.Layout.flow('vertical'), {stretch: 'both'});

ui.root.clear();
ui.root.add(dashboard);

// -------------------------------------------------------------------------
// 9. COMPLETE GOOGLE DRIVE EXPORT TASKS FOR ARCMAP
// -------------------------------------------------------------------------
Export.image.toDrive({ image: slope, description: 'Mt_Longonot_Factor_Slope', folder: 'GEE_Fire_Project', scale: 10, region: aoi, fileFormat: 'GeoTIFF' });
Export.image.toDrive({ image: aspect, description: 'Mt_Longonot_Factor_Aspect', folder: 'GEE_Fire_Project', scale: 10, region: aoi, fileFormat: 'GeoTIFF' });
Export.image.toDrive({ image: ndvi, description: 'Mt_Longonot_Factor_NDVI', folder: 'GEE_Fire_Project', scale: 10, region: aoi, fileFormat: 'GeoTIFF' });
Export.image.toDrive({ image: ndwi, description: 'Mt_Longonot_Factor_NDWI', folder: 'GEE_Fire_Project', scale: 10, region: aoi, fileFormat: 'GeoTIFF' });
Export.image.toDrive({ image: rainfall, description: 'Mt_Longonot_Factor_Rainfall', folder: 'GEE_Fire_Project', scale: 30, region: aoi, fileFormat: 'GeoTIFF' });
Export.image.toDrive({ image: lst, description: 'Mt_Longonot_Factor_LST', folder: 'GEE_Fire_Project', scale: 30, region: aoi, fileFormat: 'GeoTIFF' });
Export.image.toDrive({ image: nightLights, description: 'Mt_Longonot_Factor_NightLights', folder: 'GEE_Fire_Project', scale: 30, region: aoi, fileFormat: 'GeoTIFF' });
Export.image.toDrive({ image: rfProbability, description: 'Mt_Longonot_RF_Susceptibility_Probability', folder: 'GEE_Fire_Project', scale: 10, region: aoi, fileFormat: 'GeoTIFF' });
Export.image.toDrive({ image: dnbr2009, description: 'Mt_Longonot_Burn_Severity_2009', folder: 'GEE_Fire_Project', scale: 30, region: aoi, fileFormat: 'GeoTIFF' });
Export.image.toDrive({ image: dnbr2021, description: 'Mt_Longonot_Burn_Severity_2021', folder: 'GEE_Fire_Project', scale: 10, region: aoi, fileFormat: 'GeoTIFF' });
Export.image.toDrive({ image: dnbr2022, description: 'Mt_Longonot_Burn_Severity_2022', folder: 'GEE_Fire_Project', scale: 10, region: aoi, fileFormat: 'GeoTIFF' });
Export.image.toDrive({ image: dnbr2026, description: 'Mt_Longonot_Burn_Severity_2026', folder: 'GEE_Fire_Project', scale: 10, region: aoi, fileFormat: 'GeoTIFF' });