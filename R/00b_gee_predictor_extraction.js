// ============================================================
// Telecom Tower Siting - Ethiopia
// Predictor extraction via Google Earth Engine
// ============================================================
// Upload presence_absence.csv to your GEE Assets first
// (Assets -> New -> Table upload; name it 'tower_points').
// Set lon/lat columns when prompted. Then replace the path below.
// ============================================================

// -------------- 1. Load points and study area --------------
var points = ee.FeatureCollection('projects/YOUR-PROJECT/assets/tower_points');
var ethiopia = ee.FeatureCollection('FAO/GAUL/2015/level0')
                 .filter(ee.Filter.eq('ADM0_NAME', 'Ethiopia'));
var roi = ethiopia.geometry();

Map.centerObject(roi, 6);
Map.addLayer(ethiopia, {color: 'black'}, 'Ethiopia boundary');
Map.addLayer(points, {color: 'red'}, 'Tower points');

// -------------- 2. Elevation + terrain --------------
var dem = ee.Image('USGS/SRTMGL1_003').clip(roi);
var slope = ee.Terrain.slope(dem);
var aspect = ee.Terrain.aspect(dem);
// TRI: mean absolute difference between centre pixel and 8 neighbours
var tri = dem.subtract(dem.focalMean(1, 'square', 'pixels'))
             .abs()
             .rename('tri');

// -------------- 3. Population (WorldPop 2020 constrained) --------
var worldpop = ee.ImageCollection('WorldPop/GP/100m/pop')
                 .filter(ee.Filter.eq('country', 'ETH'))
                 .filter(ee.Filter.eq('year', 2020))
                 .first()
                 .rename('population');

// Population within 5 km (demand catchment)
var pop5km = worldpop.reduceNeighborhood({
  reducer: ee.Reducer.sum(),
  kernel: ee.Kernel.circle(5000, 'meters')
}).rename('pop_5km');

// -------------- 4. Land cover (ESA WorldCover 2021, 10 m) --------
var lc = ee.ImageCollection('ESA/WorldCover/v200').first().rename('landcover');

// -------------- 5. Built-up surface (GHSL 2020, 100 m) --------
var builtup = ee.Image('JRC/GHSL/P2023A/GHS_BUILT_S/2020').rename('built_s');

// -------------- 6. Nighttime lights (VIIRS Black Marble 2023 annual) -
var ntl = ee.ImageCollection('NOAA/VIIRS/DNB/ANNUAL_V22')
            .filter(ee.Filter.date('2023-01-01', '2023-12-31'))
            .first()
            .select('average')
            .rename('ntl');

// -------------- 7. Rainfall (CHIRPS, mean annual 2015-2023) ---------
var chirps = ee.ImageCollection('UCSB-CHG/CHIRPS/DAILY')
               .filter(ee.Filter.date('2015-01-01', '2023-12-31'))
               .sum()
               .divide(9)
               .rename('precip_mm_yr');

// -------------- 8. Water bodies (JRC Global Surface Water) ----------
var water = ee.Image('JRC/GSW1_4/GlobalSurfaceWater').select('occurrence');
var water_mask = water.gt(50).unmask(0).rename('water_mask');

// -------------- 9. Protected areas (WDPA) ---------------------------
var wdpa = ee.FeatureCollection('WCMC/WDPA/current/polygons')
             .filterBounds(roi);
var pa_raster = ee.Image(0).paint(wdpa, 1).rename('protected_area').unmask(0);

// -------------- 10. Distance layers (roads, urban, water) -----------
// Roads from OSM via MapBiomas-style proxy: we use GHSL SMOD for urban
var smod = ee.Image('JRC/GHSL/P2023A/GHS_SMOD/2020').rename('smod');

// -------------- 11. Stack all bands into one image -----------------
var stack = dem.rename('elevation')
  .addBands(slope.rename('slope'))
  .addBands(aspect.rename('aspect'))
  .addBands(tri)
  .addBands(worldpop)
  .addBands(pop5km)
  .addBands(lc)
  .addBands(builtup)
  .addBands(ntl)
  .addBands(chirps)
  .addBands(water_mask)
  .addBands(pa_raster)
  .addBands(smod);

print('Predictor stack bands:', stack.bandNames());

// -------------- 12. Sample at each point ---------------------------
var sampled = stack.sampleRegions({
  collection: points,
  properties: ['point_id', 'label', 'source', 'lon', 'lat'],
  scale: 100,     // match WorldPop resolution
  geometries: false,
  tileScale: 4    // raise if you hit memory errors
});

print('Sample preview (first 5):', sampled.limit(5));

// -------------- 13. Export to Drive as CSV -------------------------
Export.table.toDrive({
  collection: sampled,
  description: 'tower_predictors_GEE',
  fileNamePrefix: 'tower_predictors_GEE',
  fileFormat: 'CSV'
});

// -------------- 14. (Optional) Export full prediction raster -------
// Run this only when ready to map the whole country.
// Takes 30-60 min depending on tileScale.
/*
Export.image.toDrive({
  image: stack.toFloat(),
  description: 'ethiopia_predictor_stack_100m',
  scale: 100,
  region: roi,
  maxPixels: 1e13,
  crs: 'EPSG:32637'
});
*/
