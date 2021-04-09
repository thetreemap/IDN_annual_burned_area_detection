# IDN_annual_burned_area_detection
This code generates annual pre- and post-fire national (Indonesia) composite image pairs for 2019 from Sentinel-2 data in Google Earth Engine.  The data will be exported as an asset (default) to be used for the next processing (classification). To reduce the possibility of error due to the large datasets to be processed, the processing will be split into smaller grids. There are 103 grids in total (ID start from 1 to 103) covering the whole Indonesia. To process, select start and end index. Several grids can be processed at a time, but a smaller number of grids will reduce the possibility of error. These outputs can be combined again as image collection in the next step. To process other years, select the observation year.

The code can be accessed directly in GEE at the following link: https://code.earthengine.google.com/e6391f1a0dcea8c7522d2403a4d1ec41

Code created by Mohammad Agus Salim, Adrià Descals and David Gaveau from TheTreeMap (https://thetreemap.com/) 
