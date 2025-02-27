#######################################################
##### I. Citation and Abstract of the manuscript ######
#######################################################

Although textual and archaeological sources inform us about the importance of river transport in Protohistory, Antiquity, and the Middle Ages, integrating it into studies of ancient mobility remains a challenge. Empirical navigability studies are time-consuming and are only feasible for rivers well-documented by historical, archaeological, and palaeogeographical studies. This work proposes a method for realistically approximating navigable sections without empirical data by algorithmically detecting the plain sections of a river and testing its reliability as an indicator of navigability. Using 18 rivers in central-eastern Gaul, for which we have empirical knowledge of ancient navigable sections, we demonstrate that estimating the plain section of the river based on a change-point detection algorithm provides a good approximation of navigable sections. This method is applied to 48 Roman rivers where empirical information about navigable sections is scattered. A subset of these rivers is then empirically tested to validate the results obtained. Applying this method offers a new perspective on navigable areas in the Roman world, providing a reasonable first guess that could guide future empirical research into the navigability of ancient rivers.

This repository contains data and scripts used in the following paper published in JAS:

Clara Filet, C., Laroche, C., Coto-Sarmiento, M., Bongers, T. (accepted February 2025): As the water flows: a method for assessing river navigability in the past. Journal of Archaeological Science.

################################
##### II. Manuscript link ######
################################

This repository contains all the scripts for the paper "As the water flows: a method for assessing river navigability in the past" written by Clara Filet, María Coto-Sarmiento, Toon Bongers and Clément Laroche.

##########################################
##### III. Description of the files ######
##########################################

There are 5 scripts and 5 folders and each of them which consists in : 

- partsection_wo_constraint.R that applies the changepoint detection method on the elevation of points located upstream of the highest confluence location (search method M1). This file produces the partsection_wo_constraint.Rdata.
- fullsection_wo_constraint.R that applies the changepoint detection method on the elevation of all the points of the rivers (search method M2). This file produces the fullsection_wo_constraint.Rdata.
- fullsection_w_constraint.R that applies the changepoint detection method on the elevation of all the points of the rivers with the constraint that the change-points must be located upstream the confluence (search method M3). This file produces the fullsection_w_constraint.Rdata.
- Method_comp.R that produces the indicators measuring the performance on the models. It is also the script where we cluster the rivers and compare the efficiency of the methods on the different clusters. This file produces Image 11 of the paper. This file relies on the .Rdata files produced by the three first scripts. 
- Img_prod.R that produces Figures from the paper (Figures 6,7,8,9 and 10).
- the raster of the Digital Elevation Model can be found in the raster folder. 
- the navigation endpoints shapefile can be found in the empirical_endpoints folder.
- the shapefile of the riverpoints from the source to the confluence can be found in the partsection_riverpoints folder.
- the shapefile of the riverpoints from the source to the sea can be found in the fullsection_riverpoints folder. 
- all the Rdata and csv files produced by partsection_wo_constraint.R, fullsection_wo_constraint.R or fullsection_w_constraint.R are stored in the results folder.

#####################################################
##### IV. Software used and version information #####
#####################################################

R version 4.3.1 (2023-06-16 ucrt) -- "Beagle Scouts" was used to compute all results and Figures 6 to 11 of the paper. 
An API key for stadia maps is necessary to compute Figure 6. Please follow the following link explaining how to get one : https://docs.stadiamaps.com/authentication/

List of required R packages : 
- udunits2
- sf
- terra
- gridExtra
- raster
- osmdata 
- data.table
- ggplot2
- mcp
- changepoint
- segmented
- patchwork
- dplyr
- FactoMineR
- ggmap
- longitudinalData
- rPref
- reshape2

###########################
##### V. Instructions #####
###########################

Mandatory preliminary steps :
1. Download the repository as a zip.
2. Extract the repository in a folder.

To compute the results of the first search method : 
1. Open the partsection_wo_constraint.R script.
2. Run it.
3. All produced files will be found in the results folder.

To compute the results of the second search method : 
1. Open the fullsection_wo_constraint.R script.
2. Run it.
3. All produced files will be found in the results folder.

To compute the results of the third search method : 
1. Open the fullsection_w_constraint.R script.
2. Run it.
3. All produced files will be found in the results folder.

To compute the clustering of rivers according to their shape and observe the search methods results on the clusters :
1. Open the Method_comp.R script.
2. Run it.
3. All produced Figures will be found in the root of the repository.

To produce Figures 6,7,8,9,10 or/and 11 from the paper (and other interesting figures) :
1. Open the Img_prod.R script.
2. Run it.
3. All produced Figures will be found in the root of the repository.
4. A additionnal csv will be produced in the results folder. 

#######################
##### VI. Funding #####
#######################

This work was supported by the MINERVA Project (Danmarks Frie Forskningsfond (DFF) Sapere Aude research leadership grant (0163-00060B), PI: Tom Brughmans) and the Fyssen Foundation postdoctoral fellowship.

########################
##### VII. Contact #####
########################

Please, contact me if you have any questions, comments, or feedback --> mcotsar [at] gmail.com or clp.laroche[at]gmail.com

#########################
##### VIII. License #####
#########################

CC-BY 4.0