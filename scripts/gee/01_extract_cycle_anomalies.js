var samplingUnits = ee.FeatureCollection('YOUR_GEE_ASSET/sampling_unit_points') // restricted: sampling-unit coordinates must be requested from SinFlor;
Map.centerObject(samplingUnits, 7);
Map.addLayer(samplingUnits, {color: 'red'}, 'Sampling units');

var cycle1 = {name: 'cycle1', start: '2008-01-01', end: '2010-12-31'};
var cycle2 = {name: 'cycle2', start: '2014-01-01', end: '2020-12-31'};

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

function meanImage(cycle) {
  return getLST(cycle.start, cycle.end).mean().rename('lst_mean')
    .addBands(getPrecipitation(cycle.start, cycle.end).mean().rename('prec_mean'))
    .addBands(getRadiation(cycle.start, cycle.end).mean().rename('rad_mean'))
    .addBands(getNDVI(cycle.start, cycle.end).mean().rename('ndvi_mean'))
    .addBands(getRelativeHumidity(cycle.start, cycle.end).mean().rename('umid_mean'))
    .reproject({crs: 'EPSG:4326', scale: 1000});
}

var image1 = meanImage(cycle1);
var image2 = meanImage(cycle2);
var delta = image2.subtract(image1).rename([
  'delta_lst_mean',
  'delta_prec_mean',
  'delta_rad_mean',
  'delta_ndvi_mean',
  'delta_umid_mean'
]);

var sampled = delta.sampleRegions({
  collection: samplingUnits,
  scale: 1000,
  geometries: true
}).map(function(feature) {
  return feature.set('contrast', 'cycle2_minus_cycle1');
});

print('Inter-cycle environmental anomalies', sampled.limit(5));

Export.table.toDrive({
  collection: sampled,
  description: 'florestasc_intercycle_environmental_anomalies',
  fileFormat: 'CSV'
});
