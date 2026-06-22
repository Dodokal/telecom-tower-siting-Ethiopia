// ==================================================================
// Telecom Tower Siting - Ethiopia
// Predictor raster export (GeoTIFFs to Drive)
// All rasters: 100 m resolution, EPSG:32637, clipped to Ethiopia
// ==================================================================
// Run the whole script -> go to Tasks tab -> click RUN on each task.
// Expect 15-40 minutes per raster depending on size.
// Create a Google Drive folder called "ETH_towers" beforehand.
// ==================================================================

// -------------- 0. Study area ---------------------------------------
var ethiopia = ee.FeatureCollection('FAO/GAUL/2015/level0')
                 .filter(ee.Filter.eq('ADM0_NAME', 'Ethiopia'));
var roi = ethiopia.geometry();
Map.centerObject(roi, 6);
Map.addLayer(ethiopia, {color: 'black'}, 'Ethiopia');

var EXPORT_CRS   = 'EPSG:32637';   // UTM 37N (covers most of Ethiopia)
var EXPORT_SCALE = 100;            // metres
var FOLDER       = 'ETH_towers';   // create this folder in Drive first

// Helper to keep export calls tidy
function exportImg(img, name) {
  Export.image.toDrive({
    image: img.toFloat().clip(roi),
    description: name,
    folder: FOLDER,
    fileNamePrefix: name,
    region: roi,
    scale: EXPORT_SCALE,
    crs: EXPORT_CRS,
    maxPixels: 1e13
  });
}

// ==================================================================
// TERRAIN GROUP (derived from SRTM 30 m)
// ==================================================================
var dem    = ee.Image('USGS/SRTMGL1_003');
var slope  = ee.Terrain.slope(dem);
var aspect = ee.Terrain.aspect(dem);
// TRI (Riley et al.): mean absolute diff between centre pixel and its 8 neighbours
var tri = dem.subtract(dem.focalMean(1, 'square', 'pixels')).abs();
// TPI (Topographic Position Index): pixel minus mean of 500 m neighbourhood
var tpi = dem.subtract(
  dem.reduceNeighborhood({
    reducer: ee.Reducer.mean(),
    kernel: ee.Kernel.circle(500, 'meters')
  })
);

exportImg(dem,    'ETH_elevation');
exportImg(slope,  'ETH_slope');
exportImg(aspect, 'ETH_aspect');
exportImg(tri,    'ETH_tri');
exportImg(tpi,    'ETH_tpi');

// ==================================================================
// DEMAND GROUP (population, settlement, lights)
// ==================================================================

// WorldPop 2020 constrained, Ethiopia
var pop = ee.ImageCollection('WorldPop/GP/100m/pop')
            .filter(ee.Filter.eq('country', 'ETH'))
            .filter(ee.Filter.eq('year', 2020))
            .first();

// Population summed within 5 km neighbourhood (demand catchment)
var pop5km = pop.reduceNeighborhood({
  reducer: ee.Reducer.sum(),
  kernel: ee.Kernel.circle(5000, 'meters')
});

exportImg(pop,    'ETH_population');
exportImg(pop5km, 'ETH_population_5km');

// VIIRS Black Marble annual nighttime lights 2023
var ntl = ee.ImageCollection('NOAA/VIIRS/DNB/ANNUAL_V22')
            .filter(ee.Filter.date('2023-01-01', '2023-12-31'))
            .first()
            .select('average');
exportImg(ntl, 'ETH_nightlights_2023');

// GHSL Built-up surface 2020 (m^2 per cell) - urbanisation intensity
var built = ee.Image('JRC/GHSL/P2023A/GHS_BUILT_S/2020');
exportImg(built, 'ETH_builtup_2020');

// GHSL Settlement Model 2020 - degree of urbanisation (categorical)
var smod = ee.Image('JRC/GHSL/P2023A/GHS_SMOD/2020');
exportImg(smod, 'ETH_smod_2020');

// ==================================================================
// ENVIRONMENT / LAND GROUP
// ==================================================================

// ESA WorldCover 2021 - land cover (10 m native, resampled to 100 m)
var lc = ee.ImageCollection('ESA/WorldCover/v200').first();
exportImg(lc, 'ETH_landcover_2021');

// CHIRPS mean annual rainfall 2015-2023
var precip = ee.ImageCollection('UCSB-CHG/CHIRPS/DAILY')
               .filter(ee.Filter.date('2015-01-01', '2023-12-31'))
               .sum()
               .divide(9);
exportImg(precip, 'ETH_precip_annual_mm');

// JRC permanent water occurrence (%)
var water = ee.Image('JRC/GSW1_4/GlobalSurfaceWater').select('occurrence');
exportImg(water.unmask(0), 'ETH_water_occurrence');

// WDPA protected areas (binary raster: 1 = protected, 0 = not)
var wdpa = ee.FeatureCollection('WCMC/WDPA/current/polygons').filterBounds(roi);
var pa   = ee.Image(0).paint(wdpa, 1).unmask(0).rename('protected');
exportImg(pa, 'ETH_protected_areas');

// ==================================================================
// DERIVED DISTANCE-TO-WATER RASTER
// ==================================================================
// Distance from each pixel to nearest permanent water pixel (metres)
var waterBinary = water.gt(50).selfMask();
var distToWater = waterBinary.fastDistanceTransform(5000).sqrt()
                             .multiply(ee.Image.pixelArea().sqrt());
// Note: fastDistanceTransform returns squared pixel distance.
// For a clean metric version, use the cumulative cost approach below instead:
var distToWater2 = waterBinary.fastDistanceTransform().sqrt()
                              .multiply(EXPORT_SCALE)
                              .rename('dist_water_m');
exportImg(distToWater2, 'ETH_dist_to_water');

// ==================================================================
// WHAT'S NOT HERE (download outside GEE):
//   - OSM roads + distance-to-road raster  (use Geofabrik + R/QGIS)
//   - Predictive electricity grid          (Arderne et al., Zenodo 5815142)
//   - Distance to existing towers          (build from your OSM set in R)
//   - Ookla mobile performance tiles       (GitHub teamookla/ookla-open-data)
//   - Lightning density (optional)         (NASA LIS/OTD HRFC NetCDF)
// ==================================================================

print('Tasks created. Open the Tasks tab (top-right) and hit RUN on each.');
print('Tip: run terrain + population tasks first; they are the largest.');
