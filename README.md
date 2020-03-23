# eBird Analysis WoodThrush

This repository contains the code and data to reproduce the analyses in a paper currently in review: 
Johnston, Hochachka, Strimas-Mackey, Ruiz Gutierrez, Miller, Auer, Kelling, Fink. Analytical guidelines to increase the value of citizen science data: using eBird data to estimate species occurrence, (In review)

For a more general and comprehensive guide to analysing eBird data, including more code, comments, and explanations, we recommend this alternative source of code: 
Strimas-Mackey, Hochachka, Ruiz-Gutierrez, Robinson, Miller, Auer, Kelling, Fink, Johnston. 2020. Best Practices for Using eBird Data. Version 1.0. https://cornelllabofornithology.github.io/ebird-best-practices/. Cornell Lab of Ornithology, Ithaca, New York. https://doi.org/10.5281/zenodo.3620739


## Processing the data

The first few scripts contain data to process the code. The processed datasets have been saved within the `data_proc/` folder. So if you want to just run the models, you can proceed directly to scripts 05, 06, and 07.  

`00_gis-data.R` 			
Read in and prepare the BCR boundaries and mapping layers
	
`01_ebird-data.R` 		
Filter the eBird data from a local version of the eBird Basic Dataset (EBD)
select certain species, region, season. 

In order to run this, you will need to download your own local version of the EBD and use `auk_set_ebd_path` to define the folder where this is located. BUT, the full dataset requires hundreds of GB of space and once you have downloaded this it takes 3-4 hours to run `01_ebird-data.R`. So proceed with caution! 

`02_identify_bbs.R` 		
Identify the BBS routes and stops that are within the eBird dataset

`03_validation_data.R` 	
Split the data into training and validation datasets, based on year, BBS status, etc. 

`04_habitat-covariates.R`	
Extract and process the MODIS landcover information for each checklist

Some complex setup can be required for the MODIS data processing. We recommend that if you want to run this for your own data, you follow the instructions in https://cornelllabofornithology.github.io/ebird-best-practices/


## Running the models

`05_a_25_encounter_models.R`
Run 25 sets, each with 7 encounter rate random forest models. 

`05_b_encounter_models_plot.R`
Run a single set of 7 encounter rate models and create maps

`06_a_occupancy_models_plot.R` 			
Run a single set of 6 occupancy models

`07_sample_size_encounter_models.R`
Run 25 sets, each with 2 encounter rate random forest models and 5 different sample sizes


## Functions

The `R/` folder contains several functions used to help run the models. They most important/relevant functions are: 

`fit_model_enc` 
Runs a single encounter rate model, taking as parameters instructions about how to subsample the data, what type of model to run (maxent or random forest) and whether to fit covariates

`predict_model_enc` 
Takes the output from a fit_model_enc run and predicts to a given dataset

`fit_model_occu`
Runs a single occupancy model, taking as as parameters instructions about how to subsample the data and whether to fit covariates

`validate`
Takes the output from a fit_model_enc run, predicts to a new dataset, then calculates performance metrics. 


## Directory Structure

The code is setup to read most 'data' files for modelling from `data_proc`
Results will be put within subfolders into the folders `figures/` and `output/` 
`R` contains functions used in the modelling
