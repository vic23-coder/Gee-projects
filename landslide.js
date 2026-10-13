// =========================================================================
// 1. AOI DEFINITION & BOUNDARY BUFFER FIX
// =========================================================================
var aoi = ee.FeatureCollection('projects/northern-anchor-441619-g4/assets/SIGOR'); 
var aoiGeom = aoi.geometry();
var aoiGeomBuffered = aoiGeom.buffer(5000);

print('Buffered AOI area (sq km):', aoiGeomBuffered.area().divide(1e6));

// =========================================================================
// 2. SENTINEL-2 CLOUD MASK FUNCTION
// =========================================================================
function maskS2clouds(image) {
  var qa = image.select('QA60');
  var cloudBitMask = 1 << 10;
  var cirrusBitMask = 1 << 11;
  var mask = qa.bitwiseAnd(cloudBitMask).eq(0)
      .and(qa.bitwiseAnd(cirrusBitMask).eq(0));
  return image.updateMask(mask);
}

// =========================================================================
// 3. SENTINEL-2 IMAGERY (POST-EVENT)
// =========================================================================
var sentinel2020 = ee.ImageCollection('COPERNICUS/S2_SR_HARMONIZED')
  .filterBounds(aoiGeomBuffered)
  .filterDate('2020-04-18', '2020-07-30')
  .filter(ee.Filter.lt('CLOUDY_PIXEL_PERCENTAGE', 50))
  .map(maskS2clouds)
  .median()
  .clip(aoiGeomBuffered);

var visParams = {bands: ['B4', 'B3', 'B2'], min: 0, max: 3000};

// =========================================================================
// 3B. SENTINEL-2 IMAGERY (PRE-EVENT)
// =========================================================================
var sentinel2019 = ee.ImageCollection('COPERNICUS/S2_SR_HARMONIZED')
  .filterBounds(aoiGeomBuffered)
  .filterDate('2019-04-18', '2019-07-30')
  .filter(ee.Filter.lt('CLOUDY_PIXEL_PERCENTAGE', 50))
  .map(maskS2clouds)
  .median()
  .clip(aoiGeomBuffered);

Map.addLayer(sentinel2019, visParams, 'Sentinel-2 2019 (Pre-Event)');

// =========================================================================
// 4. KNOWN LANDSLIDE POINTS (CLASS = 1)
// =========================================================================
var landslidePoints = geometry.map(function(f) {
  return f.set('class', 1);
});

print('Landslide points count:', landslidePoints.size());

// =========================================================================
// 5. NON-LANDSLIDE SAMPLING (CLASS = 0)
// =========================================================================
var landslideBuffer = landslidePoints.geometry().buffer(300);
var samplingArea = aoiGeomBuffered.difference(landslideBuffer, 1);

var nonLandslidePoints = ee.FeatureCollection.randomPoints({
  region: samplingArea,
  points: landslidePoints.size(),
  seed: 42
}).map(function(f) {
  return f.set('class', 0);
});

var inventory = landslidePoints.merge(nonLandslidePoints);
print('Total Inventory Points:', inventory.size());

// =========================================================================
// 6. MAP DISPLAY & VISUALIZATION
// =========================================================================
Map.centerObject(landslidePoints, 11);
Map.addLayer(sentinel2020, visParams, 'Sentinel-2 (Cloud-Masked Composite)');
Map.addLayer(aoiGeomBuffered, {color: 'white'}, 'Sigor AOI (5km Buffer)', false);
Map.addLayer(landslideBuffer, {color: 'yellow'}, '300m Exclusion Buffer');
Map.addLayer(nonLandslidePoints, {color: 'blue'}, 'Non-Landslide Points (0)');
Map.addLayer(landslidePoints, {color: 'red'}, 'Landslide Points (1)');

// =========================================================================
// CONDITIONING FACTOR: ELEVATION & SLOPE
// =========================================================================
var dem = ee.Image('USGS/SRTMGL1_003').clip(aoiGeomBuffered);
var slope = ee.Terrain.slope(dem);

Map.addLayer(dem, {min: 900, max: 3400, palette: ['blue','green','yellow','brown','white']}, 'Elevation');
Map.addLayer(slope, {min: 0, max: 60, palette: ['green','yellow','orange','red']}, 'Slope (degrees)');

print('Elevation range:', dem.reduceRegion({
  reducer: ee.Reducer.minMax(), geometry: aoiGeomBuffered, scale: 30, maxPixels: 1e9
}));

print('Slope range:', slope.reduceRegion({
  reducer: ee.Reducer.minMax(), geometry: aoiGeomBuffered, scale: 30, maxPixels: 1e9
}));

// =========================================================================
// CONDITIONING FACTOR: RAINFALL (CHIRPS)
// =========================================================================
var chirps = ee.ImageCollection('UCSB-CHG/CHIRPS/DAILY')
  .filterDate('2013-01-01', '2023-01-01')
  .filterBounds(aoiGeomBuffered);

var annualRainfall = chirps.select('precipitation').sum()
  .divide(10)
  .clip(aoiGeomBuffered)
  .rename('rainfall');

Map.addLayer(annualRainfall, {min: 500, max: 2000, palette: ['white','blue','darkblue']}, 'Mean Annual Rainfall (mm)');

print('Rainfall range (mm/year):', annualRainfall.reduceRegion({
  reducer: ee.Reducer.minMax(), geometry: aoiGeomBuffered, scale: 5000, maxPixels: 1e9
}));

// =========================================================================
// CONDITIONING FACTOR: NDVI
// =========================================================================
var ndvi = sentinel2019.normalizedDifference(['B8', 'B4']).rename('NDVI');

Map.addLayer(ndvi, {min: -0.2, max: 0.8, palette: ['brown','yellow','green','darkgreen']}, 'NDVI (Vegetation Cover)');

print('NDVI range:', ndvi.reduceRegion({
  reducer: ee.Reducer.minMax(), geometry: aoiGeomBuffered, scale: 30, maxPixels: 1e9
}));

// =========================================================================
// CONDITIONING FACTOR: LAND USE / LAND COVER (LULC)
// =========================================================================
var worldCover = ee.ImageCollection('ESA/WorldCover/v200')
  .first()
  .clip(aoiGeomBuffered)
  .rename('lulc');

var lulcPalette = [
  '006400','ffbb22','ffff4c','f096ff','fa0000',
  'b4b4b4','f0f0f0','0064c8','0096a0','00cf75','fae6a0'
];

Map.addLayer(worldCover, {min: 10, max: 100, palette: lulcPalette}, 'Land Cover (ESA WorldCover)');

print('LULC classes present:', worldCover.reduceRegion({
  reducer: ee.Reducer.frequencyHistogram(), geometry: aoiGeomBuffered, scale: 30, maxPixels: 1e9
}));

// =========================================================================
// CONDITIONING FACTOR: SOIL TEXTURE
// =========================================================================
var soilTexture = ee.Image('OpenLandMap/SOL/SOL_TEXTURE-CLASS_USDA-TT_M/v02')
  .select('b0')
  .clip(aoiGeomBuffered)
  .rename('soil_texture');

Map.addLayer(soilTexture, {min: 1, max: 12, palette: ['8B4513','A0522D','CD853F','DEB887','F4A460','D2691E','BC8F8F','FFDEAD','FFE4B5','FFF8DC','FFFACD','F5DEB3']}, 'Soil Texture Class');

print('Soil texture classes present:', soilTexture.reduceRegion({
  reducer: ee.Reducer.frequencyHistogram(), geometry: aoiGeomBuffered, scale: 250, maxPixels: 1e9
}));

// =========================================================================
// CONDITIONING FACTOR: GROUNDWATER (NASA/GLDAS)
// =========================================================================
var gldas = ee.ImageCollection('NASA/GLDAS/V022/CLSM/G025/DA1D')
  .filterDate('2020-01-01', '2020-12-31')
  .filterBounds(aoiGeomBuffered);

var groundwater = gldas.select('GWS_tavg').mean()
  .clip(aoiGeomBuffered)
  .rename('groundwater');

Map.addLayer(
  groundwater,
  {
    min: 550,
    max: 700,
    palette: ['red','yellow','green','blue']
  },
  'Groundwater Storage'
);

print('Groundwater range:', groundwater.reduceRegion({
  reducer: ee.Reducer.minMax(), geometry: aoiGeomBuffered, scale: 27830, maxPixels: 1e9
}));

// =========================================================================
// CONDITIONING FACTOR: SOIL MOISTURE (NASA/FLDAS)
// =========================================================================
var fldas = ee.ImageCollection('NASA/FLDAS/NOAH01/C/GL/M/V001')
  .filterDate('2020-01-01', '2020-12-31')
  .filterBounds(aoiGeomBuffered);

var soilMoisture = fldas.select('SoilMoi00_10cm_tavg').mean()
  .clip(aoiGeomBuffered)
  .rename('soil_moisture');

Map.addLayer(
  soilMoisture,
  {
    min: 0.279,
    max: 0.387,
    palette: ['brown','yellow','green','blue']
  },
  'Soil Moisture'
);

print('Soil moisture range:', soilMoisture.reduceRegion({
  reducer: ee.Reducer.minMax(), geometry: aoiGeomBuffered, scale: 11132, maxPixels: 1e9
}));

// =========================================================================
// PHASE B: RECLASSIFICATION - CONVERT EACH FACTOR TO 1-5 SUSCEPTIBILITY SCALE
// (Groundwater & Soil Moisture reclass added once we see their printed ranges)
// =========================================================================

// --- 1. SLOPE (FIXED: moderate slopes = highest susceptibility, matching paper's finding) ---
var slopeReclass = ee.Image(1)
  .where(slope.gte(2).and(slope.lt(5)), 3)
  .where(slope.gte(5).and(slope.lt(15)), 5)   // <-- moderate slope = HIGHEST susceptibility now
  .where(slope.gte(15).and(slope.lt(25)), 4)
  .where(slope.gte(25).and(slope.lt(40)), 3)
  .where(slope.gte(40), 2)                     // very steep = lower (sheds material, less buildup)
  .clip(aoiGeomBuffered)
  .rename('slope_reclass');
// --- 2. ELEVATION ---
var elevReclass = ee.Image(1)
  .where(dem.gte(1200).and(dem.lt(1800)), 2)
  .where(dem.gte(1800).and(dem.lt(2400)), 3)
  .where(dem.gte(2400).and(dem.lt(3000)), 4)
  .where(dem.gte(3000), 5)
  .clip(aoiGeomBuffered)
  .rename('elev_reclass');

// --- 3. RAINFALL ---
var rainfallReclass = ee.Image(1)
  .where(annualRainfall.gte(800).and(annualRainfall.lt(1000)), 2)
  .where(annualRainfall.gte(1000).and(annualRainfall.lt(1200)), 3)
  .where(annualRainfall.gte(1200).and(annualRainfall.lt(1400)), 4)
  .where(annualRainfall.gte(1400), 5)
  .clip(aoiGeomBuffered)
  .rename('rainfall_reclass');

// --- 4. NDVI ---
var ndviReclass = ee.Image(5)
  .where(ndvi.gte(0).and(ndvi.lt(0.2)), 4)
  .where(ndvi.gte(0.2).and(ndvi.lt(0.4)), 3)
  .where(ndvi.gte(0.4).and(ndvi.lt(0.6)), 2)
  .where(ndvi.gte(0.6), 1)
  .clip(aoiGeomBuffered)
  .rename('ndvi_reclass');

// --- 5. LULC ---
var lulcReclass = ee.Image(3)
  .where(worldCover.eq(10), 1)
  .where(worldCover.eq(20), 2)
  .where(worldCover.eq(30), 3)
  .where(worldCover.eq(40), 4)
  .where(worldCover.eq(50), 5)
  .where(worldCover.eq(60), 5)
  .clip(aoiGeomBuffered)
  .rename('lulc_reclass');

// --- 6. SOIL TEXTURE ---
var soilReclass = ee.Image(3)
  .where(soilTexture.eq(1).or(soilTexture.eq(2)).or(soilTexture.eq(3)), 5)
  .where(soilTexture.eq(4).or(soilTexture.eq(5)).or(soilTexture.eq(6)), 4)
  .where(soilTexture.eq(7).or(soilTexture.eq(8)), 3)
  .where(soilTexture.eq(9).or(soilTexture.eq(10)), 2)
  .where(soilTexture.eq(11).or(soilTexture.eq(12)), 1)
  .clip(aoiGeomBuffered)
  .rename('soil_reclass');
// --- 7. GROUNDWATER (higher storage = wetter subsurface = more susceptible) ---
var groundwaterReclass = ee.Image(1)
  .where(groundwater.gte(580).and(groundwater.lt(610)), 2)
  .where(groundwater.gte(610).and(groundwater.lt(640)), 3)
  .where(groundwater.gte(640).and(groundwater.lt(670)), 4)
  .where(groundwater.gte(670), 5)
  .clip(aoiGeomBuffered)
  .rename('groundwater_reclass');


// --- 8. SOIL MOISTURE ---
// Higher soil moisture = more saturation = higher susceptibility

var soilMoistureReclass = ee.Image(1)
  .where(soilMoisture.gte(0.301).and(soilMoisture.lt(0.322)), 2)
  .where(soilMoisture.gte(0.322).and(soilMoisture.lt(0.344)), 3)
  .where(soilMoisture.gte(0.344).and(soilMoisture.lt(0.365)), 4)
  .where(soilMoisture.gte(0.365), 5)
  .clip(aoiGeomBuffered)
  .rename('soil_moisture_reclass');
// =========================================================================
// VISUALIZE ONE AS A CHECK
// =========================================================================
Map.addLayer(slopeReclass, {min: 1, max: 5, palette: ['green','yellow','orange','red','darkred']}, 'Slope Reclassified (1-5)');

// =========================================================================
// SAMPLE ALL 8 FACTORS AT EVERY INVENTORY POINT (landslide + non-landslide)
// =========================================================================

// Stack all 8 reclassified factors into one multi-band image
var allFactors = slope.rename('slope')
  .addBands(dem.rename('elevation'))
  .addBands(annualRainfall.rename('rainfall'))
  .addBands(ndvi.rename('ndvi'))
  .addBands(worldCover.rename('lulc'))
  .addBands(soilTexture.rename('soil_texture'))
  .addBands(groundwater.rename('groundwater'))
  .addBands(soilMoisture.rename('soil_moisture'));

// Sample the raw (non-reclassified) values at each of your 60 points
var sampledPoints = allFactors.sampleRegions({
  collection: inventory,
  properties: ['class'],
  scale: 30,
  geometries: true
});

print('Sampled values at all points:', sampledPoints);

// =========================================================================
// COMPARE AVERAGE VALUES: LANDSLIDE POINTS VS NON-LANDSLIDE POINTS
// =========================================================================

var landslideOnly = sampledPoints.filter(ee.Filter.eq('class', 1));
var nonLandslideOnly = sampledPoints.filter(ee.Filter.eq('class', 0));

print('AVERAGE VALUES AT LANDSLIDE POINTS:', landslideOnly.reduceColumns({
  reducer: ee.Reducer.mean().repeat(8),
  selectors: ['slope', 'elevation', 'rainfall', 'ndvi', 'lulc', 'soil_texture', 'groundwater', 'soil_moisture']
}));

print('AVERAGE VALUES AT NON-LANDSLIDE POINTS:', nonLandslideOnly.reduceColumns({
  reducer: ee.Reducer.mean().repeat(8),
  selectors: ['slope', 'elevation', 'rainfall', 'ndvi', 'lulc', 'soil_texture', 'groundwater', 'soil_moisture']
}));

// =========================================================================
// 7. CREATE THE LANDSLIDE SUSCEPTIBILITY INDEX (WEIGHTED SUM) - REVISED
// =========================================================================
// Weights adjusted based on Sigor-specific point-sampling evidence:
// Slope and rainfall showed real separation between landslide/non-landslide
// points, so both got the highest weight. LULC and soil texture showed
// weak but present separation, kept moderate. NDVI showed a REVERSED
// pattern (landslide points had MORE vegetation), so it's down-weighted
// heavily. Groundwater and soil moisture showed almost no separation
// (a known limitation from their coarse ~10-25km resolution relative to
// the Sigor AOI), so both are down-weighted to near-zero.

var susceptibilityIndex = slopeReclass.multiply(0.30)
  .add(rainfallReclass.multiply(0.30))
  .add(elevReclass.multiply(0.15))
  .add(lulcReclass.multiply(0.10))
  .add(soilReclass.multiply(0.08))
  .add(ndviReclass.multiply(0.05))
  .add(groundwaterReclass.multiply(0.01))
  .add(soilMoistureReclass.multiply(0.01))
  .clip(aoiGeomBuffered)
  .rename('susceptibility');

// Display susceptibility index continuous scale
Map.addLayer(
  susceptibilityIndex,
  {
    min: 1,
    max: 5,
    palette: ['#2b83ba', '#abdda4', '#ffffbf', '#fdae61', '#d7191c']
  },
  'Landslide Susceptibility Index (Continuous)'
);

// Compute min, max, and percentiles to inform final classification thresholds
var stats = susceptibilityIndex.reduceRegion({
  reducer: ee.Reducer.percentile([20, 40, 60, 80]),
  geometry: aoiGeom,
  scale: 30,
  maxPixels: 1e9
});
print('Susceptibility Index Quantile Thresholds (20%, 40%, 60%, 80%):', stats);

// =========================================================================
// 8. RECLASSIFY SUSCEPTIBILITY INTO 5 EVENLY DISTRIBUTED CLASSES
// =========================================================================

// Equal interval or natural breakpoint classification (1: Very Low to 5: Very High)
var susceptibilityClass = ee.Image(1)
  .where(susceptibilityIndex.gte(1.8).and(susceptibilityIndex.lt(2.4)), 2)
  .where(susceptibilityIndex.gte(2.4).and(susceptibilityIndex.lt(3.0)), 3)
  .where(susceptibilityIndex.gte(3.0).and(susceptibilityIndex.lt(3.6)), 4)
  .where(susceptibilityIndex.gte(3.6), 5)
  .clip(aoiGeomBuffered)
  .rename('susceptibility_class');

// Visualization with high-contrast color scheme
var classVis = {
  min: 1,
  max: 5,
  palette: [
    '#1a9641', // 1: Very Low (Green)
    '#a6d96a', // 2: Low (Light Green)
    '#ffffbf', // 3: Moderate (Yellow)
    '#fdae61', // 4: High (Orange)
    '#d7191c'  // 5: Very High (Red)
  ]
};

// Display final classified map
Map.addLayer(susceptibilityClass, classVis, 'Final Landslide Susceptibility Classes');

// =========================================================================
// 9. VALIDATION & ACCURACY CHECK WITH INVENTORY POINTS
// =========================================================================

// Sample the susceptibility class at your known landslide points (class = 1)
var landslideValidation = susceptibilityClass.sampleRegions({
  collection: landslidePoints,
  properties: ['class'],
  scale: 30
});

print('Landslide Points - Susceptibility Class Distribution:', 
  landslideValidation.reduceColumns({
    reducer: ee.Reducer.frequencyHistogram(),
    selectors: ['susceptibility_class']
  })
);

// =========================================================================
// 9. LEGEND
// =========================================================================

var legend = ui.Panel({
  style: {
    position: 'bottom-left',
    padding: '8px 15px'
  }
});

var legendTitle = ui.Label({
  value: 'Landslide Susceptibility',
  style: {
    fontWeight: 'bold',
    fontSize: '14px',
    margin: '0 0 6px 0'
  }
});

legend.add(legendTitle);

var classes = [
  'Very Low',
  'Low',
  'Moderate',
  'High',
  'Very High'
];

var colors = [
  'green',
  'lightgreen',
  'yellow',
  'orange',
  'red'
];

for (var i = 0; i < 5; i++) {

  var colorBox = ui.Label({
    style: {
      backgroundColor: colors[i],
      padding: '8px',
      margin: '0 0 4px 0'
    }
  });

  var description = ui.Label({
    value: classes[i],
    style: {
      margin: '0 0 4px 6px'
    }
  });

  legend.add(
    ui.Panel({
      widgets: [colorBox, description],
      layout: ui.Panel.Layout.Flow('horizontal')
    })
  );
}

Map.add(legend);

print('AVERAGE VALUES AT LANDSLIDE POINTS:', landslideOnly.reduceColumns({
  reducer: ee.Reducer.mean().repeat(8),
  selectors: ['slope', 'elevation', 'rainfall', 'ndvi', 'lulc', 'soil_texture', 'groundwater', 'soil_moisture']
}));

print('AVERAGE VALUES AT NON-LANDSLIDE POINTS:', nonLandslideOnly.reduceColumns({
  reducer: ee.Reducer.mean().repeat(8),
  selectors: ['slope', 'elevation', 'rainfall', 'ndvi', 'lulc', 'soil_texture', 'groundwater', 'soil_moisture']
}));

var nonLandslideValidation = susceptibilityClass.sampleRegions({
  collection: nonLandslidePoints,
  properties: ['class'],
  scale: 30
});

print('Non-Landslide Points - Susceptibility Class Distribution:', 
  nonLandslideValidation.reduceColumns({
    reducer: ee.Reducer.frequencyHistogram(),
    selectors: ['susceptibility_class']
  })
);

// =========================================================================
// PROPER VALIDATION: ROC / AUC CALCULATION
// =========================================================================

// Sample the CONTINUOUS susceptibility index (not the 1-5 classified version)
// at all 60 of your inventory points
var validationSample = susceptibilityIndex.sampleRegions({
  collection: inventory,
  properties: ['class'],
  scale: 30
});

// Pull the actual values down for calculation (safe here — only 60 points total)
validationSample.evaluate(function(fc) {
  var features = fc.features;

  var data = features.map(function(f) {
    return {
      value: f.properties.susceptibility,
      label: f.properties.class
    };
  });

  var totalPos = data.filter(function(d) { return d.label === 1; }).length;
  var totalNeg = data.filter(function(d) { return d.label === 0; }).length;

  var thresholds = [];
  for (var t = 1.0; t <= 5.0; t += 0.05) {
    thresholds.push(t);
  }

  var points = thresholds.map(function(t) {
    var tp = data.filter(function(d) { return d.label === 1 && d.value >= t; }).length;
    var fp = data.filter(function(d) { return d.label === 0 && d.value >= t; }).length;
    return { fpr: fp / totalNeg, tpr: tp / totalPos };
  });

  points.sort(function(a, b) { return a.fpr - b.fpr; });

  var auc = 0;
  for (var i = 1; i < points.length; i++) {
    var xDiff = points[i].fpr - points[i - 1].fpr;
    var yAvg = (points[i].tpr + points[i - 1].tpr) / 2;
    auc += xDiff * yAvg;
  }

  print('Estimated AUC:', auc);
});

// =========================================================================
// 10. EXPORT FINAL LANDSLIDE SUSCEPTIBILITY MAP
// =========================================================================

Export.image.toDrive({
  image: susceptibilityClass,
  description: 'Sigor_Landslide_Susceptibility',
  folder: 'GEE_Landslide',
  fileNamePrefix: 'Sigor_Landslide_Susceptibility',
  region: aoiGeomBuffered,
  scale: 30,
  maxPixels: 1e13
});

Export.image.toDrive({
  image: susceptibilityIndex,
  description: 'Sigor_Susceptibility_Index',
  folder: 'GEE_Landslide',
  fileNamePrefix: 'Sigor_Susceptibility_Index',
  region: aoiGeomBuffered,
  scale: 30,
  maxPixels: 1e13
});
Export.table.toDrive({
  collection: inventory,
  description: 'Sigor_Landslide_Inventory_Points',
  folder: 'GEE_Landslide',
  fileFormat: 'SHP'
});