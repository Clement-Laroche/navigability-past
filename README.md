This git repository contains all the scripts for the paper.

There are 5 scripts and 3 folders and each of them which consists in : 

- partsection_wo_constraint.R that applies the changepoint detection method on the elevation of points located upstream of the highest confluence location (search method M1). This file produces the partsection_wo_constraint.Rdata.
- fullsection_wo_constraint.R that applies the changepoint detection method on the elevation of all the points of the rivers (search method M2). This file produces the fullsection_wo_constraint.Rdata.
- fullsection_w_constraint.R that applies the changepoint detection method on the elevation of all the points of the rivers with the constraint that the change-points must be located upstream the confluence (search method M3). This file produces the fullsection_w_constraint.Rdata.
- Img_prod.R that produces Figures from the paper.
- Method_comp.R that produces the indicators measuring the performance on the models. It is also the script where we cluster the rivers and compare the efficiency of the methods on the different clusters. This file relies on the .Rdata files produced by the three first scripts.
- the raster can be found at the following link : https://www.earthdata.nasa.gov/sensors/srtm (this file is too large to be stored in Github).
- the navigation endpoints shapefile can be found in the empirical_endpoints folder.
- the shapefile of the riverpoints from the source to the confluence can be found in the partsection_riverpoints folder.
- the shapefile of the riverpoints from the source to the sea can be found in the fullsection_riverpoints folder. 
