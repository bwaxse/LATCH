# How to apply VA trained model (PheCode only)
To apply the LATCH phenotype for Long Covid, you may either 1) train your own model or 2) apply the VA trained model. Below are the steps for option 2), applying the VA trained model.

### Files needed from Github download:
* trained xgboost tree models:
	* `PheCode_12_inpat.json`
	* `PheCode_12_outpat.json`
	* `PheCode_3_allpat.json`
* codified feature names:
	* `xgb.model.port.clean.PheCode_12_inpat.Rdata`
	* `xgb.model.port.clean.PheCode_12_outpat.Rdata`
	* `xgb.model.port.clean.PheCode_3_allpat.Rdata`
* `apply_va_phecode_model.R` - template code to apply the model to your patient data

### Data needed from your patient population:
Your patient EHR feature data should be formatted in a table with the following columns: 
* `patient_num`: patient ID
* `u099.flag`: 0/1 flag to indicate presence of u099 icd code ever in EHR data
* `period12`: 1=pre u09.9 period, 0=post u09 period
* `inpat`: 0/1 flag to indicate inpatient(1) vs outpatient(0) at the time of SARS-CoV-2 infection (-7 to +14 days)
* 357 phecode features as listed in the codified features `.Rdata` objects. More info below. 

#### Curating PheCode (version 1.2) features
These are listed in format like `m_X041` or `n_X840.3` in the .Rdata objects listed above under "codified feature names". 
* The number found after "X" indicates the PheCode number. In these examples PheCode 041 and PheCode 840.3. 
* To find ICD codes grouped to each PheCode, refer to the PheWas catalog: https://phewascatalog.org/phewas/#phe12
* The `m_` and `n_` prefixes refer to counts curated as either
	1. `n_` -- the count of distinct dates a new-onset feature appears post-COVID-19 infection (0-6 months after infection)
 2. 		Example: If PheCode 008.5 appears on 3 different dates → n_X008.5 = 3
	3. `m_` -- the number of months (within 0-6 month window) in which a new-onset feature was observed at least once
 4. 		Example: After infection Oct, 2022, PheCode 008.5 in Nov, 2022 and March, → m_X008.5 = 2
	
### Using the code template
To run the `apply_va_phecode_model.r` script:
* install R packages `dplyr` and `xgboost`
* modify the paths in the script commented with "EDIT" flags to point to the correct paths and align model object names. 
* run the R code. 


