/*
This code is to generate pre and post disturbance annual composite and will be exported as an asset (default) to be used 
for the next processing (classification). 
To reduce the possibility of error due to large dataset that will be processed, the processing will be splitted into smaller grids. 
There are 103 grids in total (ID start from 1 to 103) covering the whole Indonesia. To process, please select start and end index. 
Several grids can be processed at a time and having smaller number of grids process will reduce the possibility of error.
These outputs can be combined again as image collection in the next step.
*/

var indexStartID=1; 
var indexEndID=7;
var boundary = ee.FeatureCollection("users/thetreemap/IDNBorneo_IslandBuffer3km"),
    index = ee.FeatureCollection("users/thetreemap/IDN_GridIndexRBI500K");
var indexAOI=index.filter(ee.Filter.gte('OBJECTID',indexStartID)).filter(ee.Filter.lte('OBJECTID',indexEndID));
// Define the AoI, intersection between grid and Indonesia's boundary with 3 km buffer.
var AOI =  boundary.geometry().intersection(indexAOI.geometry());

///////////////////////
// Composites using a moving window

var observationYear=2019;
var dt=2;//time interval
var cloudCover=80;
var stepMonth = 3
var shadowThresholdB11=700
 
/////////////////////////////
// Define Cloud mask function for S2 level2A

var bands_nd = ['B7', 'B11']
var bands_nbr = ['B8', 'B12']
var bands_rnbr = ['B12', 'B8']

var S2TOC_cloudMask = function(im){
  var cloudProbMask = im.select('MSK_CLDPRB').lt(60);
  var SCL = im.select('SCL')
  var SCL_mask = (SCL.eq(3).or(SCL.eq(1)).or(SCL.eq(8)).or(SCL.eq(9)).or(SCL.eq(10)).or(SCL.eq(11))).not()
  
  var shadow_mask = im.select(['B11']).lt(shadowThresholdB11)
    .or(im.select(['B4']).lt(200))
    
  var nd = im.normalizedDifference(bands_nd).rename('ND');
  var rnbr = im.normalizedDifference(['B12', 'B8']).rename('RNBR');
  var nbr = im.normalizedDifference(['B8', 'B12']).rename('NBR');
  var mask = cloudProbMask.and(shadow_mask.not()).and(SCL_mask);
  var water = im.select('SCL').remap(
    [1,2,3,4,5,6,7,8,9,10,11],
    [0,0,0,0,0,1,0,0,0,0 , 0])
  water = water.eq(1).rename('Water')
  
  var kernel = ee.Kernel.circle({radius: 20, units:'meters'});
  mask = mask
             .focal_min({kernel: kernel, iterations: 1})
  
  return im
          .addBands(nd)
          .addBands(rnbr)
          .addBands(nbr)
          .addBands(water)
          .updateMask(mask)
          .copyProperties(im);
}

var S2TOC_nbrMask = function(im){
  var nbrMask=im.select('NBR').lte(0)
  return im.updateMask(nbrMask);
  
}

var S2TOC_nbrMaskVeg = function(im){
  var nbrMask=im.select('NBR').gt(0)
  return im.updateMask(nbrMask);
  
}

///////////////////
// Define Dates
var listDates = ee.List.sequence(ee.Date.fromYMD(observationYear-1, 10, 1).millis(), ee.Date.fromYMD(observationYear+1, 2, 1).millis(), 84600000*dt)

var out_diff = ee.ImageCollection(listDates.map(function(dd){

  var targetDay = ee.Date(dd);
  //////////  
  // Create Composite Before dd
  var startDate = targetDay.advance(-stepMonth, 'month')
  var endDate = targetDay.advance(-1,'day')
  var S2 = ee.ImageCollection('COPERNICUS/S2_SR')
            .filterDate(startDate,endDate)
            .filterMetadata('CLOUD_COVERAGE_ASSESSMENT','less_than',cloudCover)
            .map(S2TOC_cloudMask)                
    
  var S2_before = S2.mean()
  var minNBR_1 = S2_before.normalizedDifference(bands_nbr)
  
  //////////  
  // Create Composite After dd
  startDate = targetDay
  endDate = targetDay.advance(stepMonth, 'month')//targetDay.advance(stepMonth, 'month')
  S2 = ee.ImageCollection('COPERNICUS/S2_SR')
            .filterDate(startDate,endDate)
            .filterMetadata('CLOUD_COVERAGE_ASSESSMENT','less_than',cloudCover)
            .map(S2TOC_cloudMask) 
  
  var S2_after = S2.filterDate(targetDay,targetDay.advance(1,'month')).mean() 
  var bare_S2_after = S2.map(S2TOC_nbrMask).qualityMosaic('RNBR')
  var veg_S2_after = S2.map(S2TOC_nbrMaskVeg).median()
  var combined_S2_after = new ee.ImageCollection([veg_S2_after.toInt32(),bare_S2_after.toInt32()]).mosaic()
    
  var minNBR_2 = S2_after.normalizedDifference(bands_nbr)
  
  /////////////////////
  // Make difference before and after
  var diff_minNBR = (minNBR_2.subtract(minNBR_1)).multiply(-1).rename('diff_nbr')
  var doy = ee.Image(targetDay.getRelative('day', 'year')).rename('doy')
  var out = S2_after.addBands(combined_S2_after).addBands(S2_before).addBands(diff_minNBR).addBands(doy.int());

  return out.set('system:time_start',targetDay.millis())
}))

var out_diff = out_diff.filterDate(observationYear+'-01-01',observationYear+'-12-31')


var pal = ["00810a","1bff00","fff700","ff0000"];
var visMax = {"min":0,"max":0.5,"palette":pal};
Map.addLayer(out_diff.select('diff_nbr').max().clip(AOI),visMax,'Maximum difference MinNBR',false)


// Make quality mosaics with difference image
var S2comp_diff = out_diff.qualityMosaic('diff_nbr').clip(AOI)
var waterMask = out_diff.select('Water').sum().lte(3)
S2comp_diff=S2comp_diff.updateMask(waterMask);
var imageVisParamS2_before = {"opacity":1,"bands":["B11_2","B8_2","B4_2"],"min":60,"max":5500,"gamma":1}
Map.addLayer(S2comp_diff,imageVisParamS2_before,'S2comp before disturbance',false)
var imageVisParamS2_after = {"opacity":1,"bands":["B11_1","B8_1","B4_1"],"min":60,"max":5500,"gamma":1}
Map.addLayer(S2comp_diff,imageVisParamS2_after,'S2comp after disturbance',true)


Export.image.toAsset({
  image: S2comp_diff.select(['B4_1','B5_1','B6_1','B7_1','B8_1','B8A_1','B11_1','B12_1','NBR_1']).rename(['B4','B5','B6','B7','B8','B8A','B11','B12','NBR']),
  description: 'IDNSentinel2Post'+observationYear+'-'+indexStartID+'to'+indexEndID,
  assetId: 'IDN_S2PostBurnImage'+observationYear+'/IDNSentinel2Post'+observationYear+'-'+indexStartID+'to'+indexEndID,
  scale: 20,
  region: AOI,
  maxPixels: 30000000000,
});

Export.image.toAsset({
  image: S2comp_diff.select(['B4_2','B5_2','B6_2','B7_2','B8_2','B8A_2','B11_2','B12_2','NBR_2']).rename(['B4','B5','B6','B7','B8','B8A','B11','B12','NBR']),
  description: 'IDNSentinel2Pre'+observationYear+'-'+indexStartID+'to'+indexEndID,
  assetId: 'IDN_S2PreBurnImage'+observationYear+'/IDNSentinel2Pre'+observationYear+'-'+indexStartID+'to'+indexEndID,
  scale: 20,
  region: AOI,
  maxPixels: 30000000000,
});

Map.addLayer(index,{},'Grid index',false);
Map.addLayer(AOI,{},'AOI',false);
Map.centerObject(AOI,6);
