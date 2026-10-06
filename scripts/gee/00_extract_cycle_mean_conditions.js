var samplingUnits = ee.FeatureCollection('YOUR_GEE_ASSET/sampling_unit_points') // restricted: sampling-unit coordinates must be requested from SinFlor;
Map.centerObject(samplingUnits, 7);
Map.addLayer(samplingUnits, {color: 'blue'}, 'Sampling units');

var cycles = [
  {name: 'cycle1', start: '2008-01-01', end: '2010-12-31'},
  {name: 'cycle2', start: '2014-01-01', end: '2020-12-31'}
];

function getLST(start, end) {
  return ee.ImageCollection('MODIS/061/MOD11A2')
    .filterDate(start, end)
    .select('LST_Day_1km')
    .map(function(img) {
      return img.multiply(0.02).subtract(273.15).rename('lst')
        .copyProperties(img, ['system:time_start']);
    });
}

function getPrecipitation(start, end) {
  return ee.ImageCollection('UCSB-CHG/CHIRPS/DAILY')
    .filterDate(start, end)
    .select('precipitation')
    .map(function(img) {
      return img.rename('prec').copyProperties(img, ['system:time_start']);
    });
}

function getNDVI(start, end) {
  return ee.ImageCollection('MODIS/061/MOD13Q1')
    .filterDate(start, end)
    .select('NDVI')
    .map(function(img) {
      return img.multiply(0.0001).rename('ndvi')
        .copyProperties(img, ['system:time_start']);
    });
}

function getRadiation(start, end) {
  return ee.ImageCollection('ECMWF/ERA5_LAND/DAILY_AGGR')
    .filterDate(start, end)
    .select('surface_net_solar_radiation_sum')
    .map(function(img) {
      return img.multiply(0.000001).rename('rad')
        .copyProperties(img, ['system:time_start']);
    });
}

function getRelativeHumidity(start, end) {
  var t = ee.ImageCollection('ECMWF/ERA5_LAND/DAILY_AGGR')
    .filterDate(start, end)
    .select('temperature_2m')
    .map(function(img) {
      return img.subtract(273.15).rename('t')
        .copyProperties(img, ['system:time_start']);
    });

  var td = ee.ImageCollection('ECMWF/ERA5_LAND/DAILY_AGGR')
    .filterDate(start, end)
    .select('dewpoint_temperature_2m')
    .map(function(img) {
      return img.subtract(273.15).rename('td')
        .copyProperties(img, ['system:time_start']);
    });

  var joined = ee.Join.inner().apply({
    primary: t,
    secondary: td,
    condition: ee.Filter.equals({leftField: 'system:time_start', rightField: 'system:time_start'})
  });

  return ee.ImageCollection(joined.map(function(pair) {
    var temp = ee.Image(pair.get('primary'));
    var dew = ee.Image(pair.get('secondary'));
    return dew.expression(
      '100 * (exp((17.625 * TD) / (243.04 + TD)) / exp((17.625 * T) / (243.04 + T)))',
      {TD: dew.select('td'), T: temp.select('t')}
    ).rename('umid').copyProperties(temp, ['system:time_start']);
  }));
}

function summarizeCollection(collection, prefix, includeMinMax) {
  var img = collection.mean().rename(prefix + '_mean')
    .addBands(collection.reduce(ee.Reducer.stdDev()).rename(prefix + '_sd'));
  if (includeMinMax) {
    img = img.addBands(collection.min().rename(prefix + '_min'))
      .addBands(collection.max().rename(prefix + '_max'));
  }
  return img;
}

var elevation = ee.Image('USGS/SRTMGL1_003')
  .select('elevation')
  .rename('altitude')
  .resample('bilinear')
  .reproject({crs: 'EPSG:4326', scale: 1000});

var output = ee.FeatureCollection([]);

cycles.forEach(function(cycle) {
  var lst = summarizeCollection(getLST(cycle.start, cycle.end), 'lst', true);
  var prec = summarizeCollection(getPrecipitation(cycle.start, cycle.end), 'prec', true);
  var ndvi = summarizeCollection(getNDVI(cycle.start, cycle.end), 'ndvi', true);
  var rad = summarizeCollection(getRadiation(cycle.start, cycle.end), 'rad', false);
  var umid = summarizeCollection(getRelativeHumidity(cycle.start, cycle.end), 'umid', false);

  var image = lst.addBands(prec).addBands(ndvi).addBands(rad).addBands(umid).addBands(elevation)
    .reproject({crs: 'EPSG:4326', scale: 1000});

  var sampled = image.sampleRegions({
    collection: samplingUnits,
    scale: 1000,
    geometries: true
  }).map(function(feature) {
    return feature.set('cycle', cycle.name);
  });

  output = output.merge(sampled);
});

print('Cycle-level environmental summaries', output.limit(5));

Export.table.toDrive({
  collection: output,
  description: 'florestasc_cycle_mean_environmental_conditions',
  fileFormat: 'CSV'
});
