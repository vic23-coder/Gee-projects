// ==========================================================
// COMPREHENSIVE MWEA RICE MONITORING & YIELD ESTIMATION SYSTEM
// ==========================================================

// 1. SETUP
Map.centerObject(roi, 13);
Map.setOptions('HYBRID');

// ==========================================================
// 2. FUNCTIONS: CLOUD MASKING & INDICES
// ==========================================================

function maskS2clouds(img) {
  var qa = img.select('QA60');
  var cloudMask = qa.bitwiseAnd(1 << 10).eq(0)
    .and(qa.bitwiseAnd(1 << 11).eq(0));
  return img.updateMask(cloudMask)
    .divide(10000)
    .copyProperties(img, ['system:time_start']);
}

function addComprehensiveIndices(img) {
  var ndvi = img.normalizedDifference(['B8', 'B4']).rename('NDVI');
  var ndre = img.normalizedDifference(['B8', 'B5']).rename('NDRE');
  
  var evi = img.expression(
    '2.5*((NIR-RED)/(NIR+6*RED-7.5*BLUE+1))',
    {'NIR': img.select('B8'), 'RED': img.select('B4'), 'BLUE': img.select('B2')}
  ).rename('EVI');
  
  var lswi = img.normalizedDifference(['B8', 'B11']).rename('LSWI');
  
  var savi = img.expression(
    '((NIR-RED)/(NIR+RED+0.5))*1.5',
    {'NIR': img.select('B8'), 'RED': img.select('B4')}
  ).rename('SAVI');
  
  var ndmi = img.normalizedDifference(['B8', 'B11']).rename('NDMI');
  
  var gcc = img.expression(
    'GREEN/(RED+GREEN+BLUE)',
    {'RED': img.select('B4'), 'GREEN': img.select('B3'), 'BLUE': img.select('B2')}
  ).rename('GCC');
  
  return img.addBands([ndvi, ndre, evi, lswi, savi, ndmi, gcc]);
}

// ==========================================================
// 3. DATA COLLECTION
// ==========================================================

var col = ee.ImageCollection('COPERNICUS/S2_SR_HARMONIZED')
  .filterBounds(roi)
  .filterDate('2025-06-01', '2025-11-30')
  .filter(ee.Filter.lt('CLOUDY_PIXEL_PERCENTAGE', 60))
  .map(maskS2clouds)
  .map(addComprehensiveIndices);

// Growth stages
var stage1_transplant = col.filterDate('2025-07-01', '2025-07-31').median().clip(roi);
var stage2_tillering  = col.filterDate('2025-08-15', '2025-09-15').median().clip(roi);
var stage4_heading    = col.filterDate('2025-10-15', '2025-11-15').median().clip(roi);

// ==========================================================
// 4. MASKING (BUILT-UP & RICE)
// ==========================================================

// Built-up (Dynamic World)
var dw = ee.ImageCollection('GOOGLE/DYNAMICWORLD/V1')
  .filterBounds(roi)
  .filterDate('2025-06-01', '2025-11-30');

var dwComposite = dw.select('label').reduce(ee.Reducer.mode()).clip(roi);
var builtUpMask = dwComposite.eq(6);
var nonBuiltMask = builtUpMask.eq(0);

// Rice detection
var floodingCondition = stage1_transplant.select('LSWI')
  .add(0.05)
  .gte(stage1_transplant.select('EVI'));

var matureRice = stage4_heading.select('NDVI').gte(0.6);

var riceMask = floodingCondition
  .and(matureRice)
  .and(nonBuiltMask);

// ==========================================================
// 5. YIELD ESTIMATION (NDRE-BASED)
// ==========================================================

// Empirical rice yield model
// Yield (t/ha) = 8.5 × NDRE + 1.2
var riceNDRE = stage4_heading
  .select('NDRE')
  .updateMask(riceMask);

var yieldMap = riceNDRE
  .multiply(8.5)
  .add(1.2)
  .rename('Rice_Yield_t_ha');

// Mean yield
var meanYield = yieldMap.reduceRegion({
  reducer: ee.Reducer.mean(),
  geometry: roi,
  scale: 20,
  maxPixels: 1e13
});

print('Average Rice Yield (t/ha):', meanYield);

// ==========================================================
// 6. VISUALIZATION HELPERS
// ==========================================================

var vegPalette = ['#FFD700', '#ADFF2F', '#32CD32', '#008000', '#004d00'];
var visNDVI = {min: 0.2, max: 0.9, palette: vegPalette};
var visNDRE = {min: 0.3, max: 0.65, palette: vegPalette};
var visYield = {
  min: 2,
  max: 8,
  palette: ['#fdae61', '#a6d96a', '#1a9850']
};

function overlayBuilt(img, vis) {
  var base = img.visualize(vis);
  var built = builtUpMask.visualize({palette: ['ff0000'], opacity: 1});
  return ee.ImageCollection([base, built.updateMask(builtUpMask)]).mosaic();
}

// ==========================================================
// 7. UI PANEL & LAYER SELECTOR
// ==========================================================

var panel = ui.Panel({
  style: {width: '350px', padding: '10px', position: 'bottom-right'}
});

var legendPanel = ui.Panel({
  style: {padding: '8px', border: '1px solid #ddd', margin: '10px 0'}
});

function addLegend(title, items) {
  legendPanel.clear();
  legendPanel.add(ui.Label({value: title, style: {fontWeight: 'bold'}}));
  items.forEach(function(i) {
    legendPanel.add(ui.Panel({
      widgets: [
        ui.Label({style: {backgroundColor: i.color, padding: '8px'}}),
        ui.Label(i.label)
      ],
      layout: ui.Panel.Layout.Flow('horizontal')
    }));
  });
}

var layerSelect = ui.Select({
  items: [
    {label: 'Crop Vigor (NDVI)', value: 'ndvi'},
    {label: 'Nitrogen Status (NDRE)', value: 'ndre'},
    {label: 'Estimated Rice Yield (t/ha)', value: 'yield'}
  ],
  onChange: function(val) {
    Map.layers().reset();
    
    if (val === 'ndvi') {
      Map.addLayer(overlayBuilt(stage2_tillering.select('NDVI'), visNDVI), {}, 'NDVI');
      addLegend('NDVI', [
        {color:'#FFD700', label:'Low Vigor'},
        {color:'#32CD32', label:'Healthy'},
        {color:'#004d00', label:'High Biomass'},
        {color:'ff0000', label:'Built-up'}
      ]);
    }
    
    if (val === 'ndre') {
      Map.addLayer(overlayBuilt(stage4_heading.select('NDRE'), visNDRE), {}, 'NDRE');
      addLegend('NDRE', [
        {color:'#FFD700', label:'Low Nitrogen'},
        {color:'#32CD32', label:'Adequate'},
        {color:'#004d00', label:'Excellent'},
        {color:'ff0000', label:'Built-up'}
      ]);
    }
    
    if (val === 'yield') {
      Map.addLayer(overlayBuilt(yieldMap, visYield), {}, 'Rice Yield');
      addLegend('Rice Yield (t/ha)', [
        {color:'#fdae61', label:'Low Yield'},
        {color:'#a6d96a', label:'Medium Yield'},
        {color:'#1a9850', label:'High Yield'},
        {color:'ff0000', label:'Built-up'}
      ]);
    }
  }
});

// ==========================================================
// 8. INITIALIZE
// ==========================================================

panel.add(ui.Label({
  value: 'Mwea Rice Monitoring & Yield',
  style: {fontSize: '18px', fontWeight: 'bold'}
}));

panel.add(layerSelect);
panel.add(legendPanel);

Map.add(panel);
layerSelect.setValue('yield');
