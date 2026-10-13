// =====================================================
// FLOOD MAPPING — TANA RIVER COUNTY, KENYA
// Project: Flood Mapping & Susceptibility Analysis
// Event: El Niño Floods — Peak: Nov 15–Dec 10 2023
// Methodology: Tavus et al. (2020) ISPRS Congress
// =====================================================

var galole = table;
Map.centerObject(galole, 10);
Map.addLayer(galole, {color: 'red'}, 'Tana River Boundary');

// =====================================================
// STEP 1 — CLOUD + SHADOW MASK FUNCTION
// FIX: Added cloud SHADOW masking using SCL band
// Previous version only masked clouds, not shadows
// SCL (Scene Classification Layer) in S2 SR data:
//   SCL value 3  = cloud shadow
//   SCL value 8  = medium probability cloud
//   SCL value 9  = high probability cloud
//   SCL value 10 = thin cirrus
// Combined with QA60 for maximum cloud removal
// =====================================================
function maskS2cloudsStrict(image) {
  // QA60 mask (clouds + cirrus)
  var qa = image.select('QA60');
  var cloudBitMask  = 1 << 10;
  var cirrusBitMask = 1 << 11;
  var qa_mask = qa.bitwiseAnd(cloudBitMask).eq(0)
      .and(qa.bitwiseAnd(cirrusBitMask).eq(0));

  // SCL mask (cloud shadows + clouds + cirrus)
  var scl = image.select('SCL');
  var scl_mask = scl.neq(3)   // not cloud shadow
      .and(scl.neq(8))         // not medium cloud
      .and(scl.neq(9))         // not high cloud
      .and(scl.neq(10));       // not thin cirrus

  // Apply both masks together
  return image.updateMask(qa_mask.and(scl_mask))
      .divide(10000)
      .copyProperties(image, ['system:time_start']);
}

// =====================================================
// STEP 2 — DS4: NON-FLOODED SENTINEL-2 (BASELINE)
// FIX 1: Extended window to Jan–Sep 2023
//         More scenes = better cloud-free composite
// FIX 2: Dropped cloud threshold to 10%
//         Only uses nearly cloud-free scenes
// FIX 3: Used percentile(10) instead of median()
//         10th percentile picks darkest clear pixels
//         which better represents true surface
//         reflectance without residual cloud haze
// Purpose: NDVI + MNDWI computation ONLY
// =====================================================
var DS4_collection = ee.ImageCollection('COPERNICUS/S2_SR_HARMONIZED')
  .filterBounds(galole)
  .filterDate('2023-01-01', '2023-09-30')  // wider window
  .filter(ee.Filter.lt('CLOUDY_PIXEL_PERCENTAGE', 10)) // strict
  .map(maskS2cloudsStrict);

print('DS4 Non-Flooded S2 scene count:', DS4_collection.size());

// percentile(10) = picks near-clearest pixels
// avoids bright cloud residuals that median keeps
var DS4 = DS4_collection
  .reduce(ee.Reducer.percentile([10]))
  .rename(DS4_collection.first().bandNames())
  .clip(galole);

Map.addLayer(DS4,
  {bands: ['B4', 'B3', 'B2'], min: 0, max: 0.3, gamma: 1.4},
  'DS4 — Non-Flooded RGB (Cloud-Free Baseline)');

// =====================================================
// STEP 3 — DS3: FLOODED SENTINEL-2
// Visual reference only — NOT used for any calculation
// Shown so you can see what the flood looked like
// where Sentinel-2 managed to see through clouds
// =====================================================
var DS3_collection = ee.ImageCollection('COPERNICUS/S2_SR_HARMONIZED')
  .filterBounds(galole)
  .filterDate('2023-11-15', '2023-12-10')
  .filter(ee.Filter.lt('CLOUDY_PIXEL_PERCENTAGE', 90))
  .map(maskS2cloudsStrict);

print('DS3 Flooded S2 scene count (visual only):',
  DS3_collection.size());

var DS3 = DS3_collection.mosaic().clip(galole);

Map.addLayer(DS3,
  {bands: ['B4', 'B3', 'B2'], min: 0, max: 0.3, gamma: 1.4},
  'DS3 — Flooded RGB (Visual Reference ONLY)');

// =====================================================
// STEP 4 — NDVI FROM DS4 (PAPER EQUATION 1)
// NDVI = (NIR - Red) / (NIR + Red)
// B8 = NIR (10m), B4 = Red (10m)
// Context layer: shows vegetation distribution
// Input for flood susceptibility modelling
// =====================================================
var ndvi = DS4
  .normalizedDifference(['B8', 'B4'])
  .rename('NDVI');

Map.addLayer(ndvi,
  {min: -0.2, max: 0.8,
   palette: ['brown', 'yellow', 'lightgreen', 'darkgreen']},
  'NDVI — DS4 (brown=bare/water, green=vegetation)');

// =====================================================
// STEP 5 — MNDWI FROM DS4 (PAPER EQUATION 2)
// MNDWI = (Green - SWIR) / (Green + SWIR) [Xu 2006]
// B3 = Green (10m), B11 = SWIR1 (20m)
// Shows permanent water bodies in dry season
// Used as context/feature layer — NOT flood detection
// =====================================================
var mndwi = DS4
  .normalizedDifference(['B3', 'B11'])
  .rename('MNDWI');

Map.addLayer(mndwi,
  {min: -0.5, max: 0.5,
   palette: ['8B4513', 'F5DEB3', 'FFFFFF', '40E0D0', '00008B']},
  'MNDWI — DS4 (blue=permanent water, brown=dry land)');

// =====================================================
// STEP 6 — DS1: PRE-FLOOD SENTINEL-1 SAR
// Period: Jan–Sep 2023 (extended dry season)
// FIX: Extended window to match DS4 for consistency
// IW mode, dual pol VV + VH
// C11 = VV (open water), C22 = VH (vegetation/crops)
// =====================================================
var DS1_collection = ee.ImageCollection('COPERNICUS/S1_GRD')
  .filterBounds(galole)
  .filterDate('2023-01-01', '2023-09-30')
  .filter(ee.Filter.eq('instrumentMode', 'IW'))
  .filter(ee.Filter.listContains('transmitterReceiverPolarisation', 'VV'))
  .filter(ee.Filter.listContains('transmitterReceiverPolarisation', 'VH'))
  .select(['VV', 'VH']);

print('DS1 Pre-flood SAR scene count:', DS1_collection.size());

var DS1 = DS1_collection.mean().clip(galole);

// =====================================================
// STEP 7 — DS2: PEAK FLOOD SENTINEL-1 SAR
// Period: Nov 15–Dec 10 2023 (El Niño flood peak)
// SAR sees through clouds — no cloud problem here
// This is the correct sensor for flood detection
// during El Niño persistent heavy cloud cover
// =====================================================
var DS2_collection = ee.ImageCollection('COPERNICUS/S1_GRD')
  .filterBounds(galole)
  .filterDate('2023-11-15', '2023-12-10')
  .filter(ee.Filter.eq('instrumentMode', 'IW'))
  .filter(ee.Filter.listContains('transmitterReceiverPolarisation', 'VV'))
  .filter(ee.Filter.listContains('transmitterReceiverPolarisation', 'VH'))
  .select(['VV', 'VH']);

print('DS2 Peak Flood SAR scene count:', DS2_collection.size());

var DS2 = DS2_collection.mean().clip(galole);

// =====================================================
// STEP 8 — SAR SPECKLE FILTERING
// Paper Section 3.1: Improved Lee Sigma filter 7x7
// GEE equivalent: focal_mean 70m radius
// (7 pixels × 10m = 70m)
// Applied to VV (C11) and VH (C22) separately
// Removes salt-and-pepper noise from SAR images
// =====================================================
var C11_pre  = DS1.select('VV')
  .focal_mean(70, 'circle', 'meters').rename('C11_pre');
var C11_post = DS2.select('VV')
  .focal_mean(70, 'circle', 'meters').rename('C11_post');
var C22_pre  = DS1.select('VH')
  .focal_mean(70, 'circle', 'meters').rename('C22_pre');
var C22_post = DS2.select('VH')
  .focal_mean(70, 'circle', 'meters').rename('C22_post');

// In SAR images: dark = smooth (water), bright = rough (land)
Map.addLayer(C11_pre,
  {min: -25, max: 0, palette: ['black', 'gray', 'white']},
  'C11 VV — DS1 Pre-Flood SAR (dark=water)');
Map.addLayer(C11_post,
  {min: -25, max: 0, palette: ['black', 'gray', 'white']},
  'C11 VV — DS2 Peak Flood SAR (dark=water)');
Map.addLayer(C22_pre,
  {min: -25, max: 0, palette: ['black', 'gray', 'white']},
  'C22 VH — DS1 Pre-Flood SAR');
Map.addLayer(C22_post,
  {min: -25, max: 0, palette: ['black', 'gray', 'white']},
  'C22 VH — DS2 Peak Flood SAR');

// =====================================================
// STEP 9 — SAR BACKSCATTER DIFFERENCE
// Pre-flood minus post-flood backscatter
// Flooded pixel = large DROP in backscatter
// because flat water bounces radar signal away
// Red = large drop = flooded
// White = no change = not flooded
// =====================================================
var diff_VV = C11_pre.subtract(C11_post).rename('diff_VV');
var diff_VH = C22_pre.subtract(C22_post).rename('diff_VH');

Map.addLayer(diff_VV,
  {min: -5, max: 10,
   palette: ['white', 'lightyellow', 'yellow', 'orange', 'red']},
  'VV Backscatter Difference (red=flooded)');
Map.addLayer(diff_VH,
  {min: -5, max: 10,
   palette: ['white', 'lightyellow', 'yellow', 'orange', 'red']},
  'VH Backscatter Difference (red=flooded)');

// =====================================================
// STEP 10 — SAR FLOOD MASK
// Threshold: >2dB drop = flooded
// Tana floodplain is flat so 2dB not 3dB
// VV catches open water flooding
// VH catches flooded crops and vegetation
// OR logic ensures both are captured
// focal_mode removes remaining speckle noise
// =====================================================
var sar_flood_VV = diff_VV.gt(2);
var sar_flood_VH = diff_VH.gt(2);

var sar_flood = sar_flood_VV.or(sar_flood_VH)
  .focal_mode(30, 'circle', 'meters')
  .rename('SAR_flood');

Map.addLayer(sar_flood.selfMask(),
  {palette: ['0000FF']},
  'SAR Flood Extent — PRIMARY (blue=flooded)');

// =====================================================
// STEP 11 — FINAL FLOOD vs NON-FLOOD MAP
// SAR-based only — no optical flood detection
// 1 = Flooded (red)
// 2 = Non-Flooded (green)
// Primary output for your project report
// =====================================================
var flood_map = ee.Image(2)
  .where(sar_flood.eq(1), 1)
  .clip(galole)
  .rename('flood_class');

Map.addLayer(flood_map,
  {min: 1, max: 2,
   palette: ['FF0000', '90EE90']},
  'Flood vs Non-Flood — FINAL MAP');

// =====================================================
// STEP 12 — DEM + TERRAIN
// Adjusted for flat Tana River floodplain
// Elevation: 0–150m, Slope: 0–5 degrees
// Inputs for susceptibility modelling
// =====================================================
var dem       = ee.Image('USGS/SRTMGL1_003').clip(galole);
var slope     = ee.Terrain.slope(dem);
var aspect    = ee.Terrain.aspect(dem);
var hillshade = ee.Terrain.hillshade(dem, 315, 45);

Map.addLayer(hillshade,
  {min: 150, max: 255},
  'Hillshade (enhanced for flat terrain)');
Map.addLayer(dem,
  {min: 0, max: 150,
   palette: ['006633','39B54A','E5FFCC','FFFF00','FF6600','FF0000']},
  'DEM Elevation (0–150m)');
Map.addLayer(slope,
  {min: 0, max: 5,
   palette: ['white', 'yellow', 'orange', 'red']},
  'Slope (0–5 degrees)');
Map.addLayer(aspect,
  {min: 0, max: 360,
   palette: ['blue', 'green', 'yellow', 'red', 'blue']},
  'Aspect');

// =====================================================
// STEP 13 — FLOOD AREA STATISTICS
// =====================================================
var pixelArea = ee.Image.pixelArea();
var toKm2     = 1e6;

var sar_area = sar_flood.multiply(pixelArea)
  .reduceRegion({
    reducer: ee.Reducer.sum(),
    geometry: galole.geometry(),
    scale: 10, maxPixels: 1e13
  });

var total_area = ee.Image(1).clip(galole)
  .multiply(pixelArea)
  .reduceRegion({
    reducer: ee.Reducer.sum(),
    geometry: galole.geometry(),
    scale: 10, maxPixels: 1e13
  });

print('=== FLOOD AREA STATISTICS ===');
print('Total Study Area (km²):',
  ee.Number(total_area.get('constant')).divide(toKm2));
print('Flooded Area (km²):',
  ee.Number(sar_area.get('SAR_flood')).divide(toKm2));
print('Non-Flooded Area (km²):',
  ee.Number(total_area.get('constant'))
    .subtract(ee.Number(sar_area.get('SAR_flood')))
    .divide(toKm2));

// =====================================================
// STEP 14 — EXPORTS
// Uncomment to export to Google Drive
// =====================================================

// Export.image.toDrive({
//   image: flood_map,
//   description: 'Tana_ElNino_Flood_vs_NonFlood',
//   folder: 'Flood_Tana_River',
//   region: galole.geometry(),
//   scale: 10, fileFormat: 'GeoTIFF', maxPixels: 1e13
// });

// Export.image.toDrive({
//   image: sar_flood,
//   description: 'Tana_SAR_FloodExtent',
//   folder: 'Flood_Tana_River',
//   region: galole.geometry(),
//   scale: 10, fileFormat: 'GeoTIFF', maxPixels: 1e13
// });

// Export.image.toDrive({
//   image: ndvi,
//   description: 'Tana_NDVI_Baseline',
//   folder: 'Flood_Tana_River',
//   region: galole.geometry(),
//   scale: 10, fileFormat: 'GeoTIFF', maxPixels: 1e13
// });

// Export.image.toDrive({
//   image: mndwi,
//   description: 'Tana_MNDWI_Baseline',
//   folder: 'Flood_Tana_River',
//   region: galole.geometry(),
//   scale: 10, fileFormat: 'GeoTIFF', maxPixels: 1e13
// });

// Export.image.toDrive({
//   image: dem,
//   description: 'Tana_DEM',
//   folder: 'Flood_Tana_River',
//   region: galole.geometry(),
//   scale: 30, fileFormat: 'GeoTIFF', maxPixels: 1e13
// });

// Export.image.toDrive({
//   image: slope,
//   description: 'Tana_Slope',
//   folder: 'Flood_Tana_River',
//   region: galole.geometry(),
//   scale: 30, fileFormat: 'GeoTIFF', maxPixels: 1e13
// });

// Export.image.toDrive({
//   image: aspect,
//   description: 'Tana_Aspect',
//   folder: 'Flood_Tana_River',
//   region: galole.geometry(),
//   scale: 30, fileFormat: 'GeoTIFF', maxPixels: 1e13
// });