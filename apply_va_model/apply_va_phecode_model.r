library(dplyr)
library(xgboost)
rm(list=ls())

### import pt data
load(file="my_path/patient.data.Rdata") ### EDIT TO PATH TO PATIENT DATA OBJECT. referenced as "my_patient_data" variable below.

#### import the xgboost jsons ####
indir = "U:/My Documents/Phenotyping/PASC/xgboost_models" ### EDIT TO DIRECTORY CONTAINING XGBOOST MODEL FILES
fnames = c("PheCode_12_inpat.json", "PheCode_12_outpat.json", "PheCode_3_allpat.json")
xgb.model = list()
for (i in 1:length(fnames)) {
    f = file.path(indir, fnames[i])
    m_name = sub(".json", "", fnames[i])
    m = xgb.load(f)
    xgb.model[[m_name]] = m
}

### import codified feature names ###
indir = "U:/My Documents/Phenotyping/PASC/feature_names" ### EDIT TO DIRECTORY CONTAINING FEATURE NAME FILES
fnames = c("xgb.model.port.clean.PheCode_12_inpat.Rdata", 
           "xgb.model.port.clean.PheCode_12_outpat.Rdata", 
           "xgb.model.port.clean.PheCode_3_allpat.Rdata")
feature.names = list()
for (i in 1:length(fnames)) {
  f = file.path(indir, fnames[i])
  m_name = names(xgb.model)[i]
  load(f)
  feature.names[[m_name]] = xx$feature_names
}

### apply xgboost models to pt data to get xgboost scores
xgb.pred0=NULL
for (i in 1:length(xgb.model)){
    m = xgb.model[[i]]
    m_name = names(xgb.model)[[i]]
    print(m_name)

    ### get correct feature list
    feature.sel.keep <- feature.names[[i]] 
    if(grepl("allpat", m_name)){
        feature.sel.keep=c("inpat", feature.sel.keep)
    }
    dtest = xgb.DMatrix(
      data=data.matrix(my_patient_data[,which(colnames(my_patient_data) %in% feature.sel.keep)]),
      label=my_patient_data[,"u099.flag"]
    )
    xgb.pred0=cbind(xgb.pred0, predict(m, dtest))
}
colnames(xgb.pred0) = ls(xgb.model)


### clean up dataset for application of regression model

### variables needed:
## - patient_num: patient ID
## - u099.flag: 0/1 flag to indicate presence of u099 icd code
## - period12: 1=pre u09.9 period, 0=post u09 period.
## - inpat: 0/1 flag to indicate inpatient(1) vs outpatient(0)
xgb.pred = data.frame(patient_num = my_patient_data$patient_num,
                      u099.flag = my_patient_data$u099.flag,
                      period12 = my_patient_data$period12,
                      inpat = my_patient_data$inpat,
                      xgb.pred0)


## "cohort alignment step": align the xgboost scores based on inpat/outpat status, infection period
##  aligned xgboost score is called "model.score"
dat.ssl <- data.frame(
  U099_Count = my_patient_data$U099_Count,
  xgb.pred
)
dat.ssl <- dat.ssl %>% 
  mutate(model.score = case_when(
    period12==1 & inpat==1 ~ PheCode_12_inpat,
    period12==1 & inpat==0 ~ PheCode_12_outpat,
    period12==0 ~ PheCode_3_allpat
  ))

## Application of regression model
dat.ssl$ssl_ord=g.logit(0.6372599*log(dat.ssl$U099_Count+1)+0.2568035*dat.ssl$model.score)
dat.ssl$ssl_who1=g.logit(1.4151418*log(dat.ssl$U099_Count+1)+0.301594*dat.ssl$model.score)
dat.ssl$ssl_who=g.logit(1.0843817*log(dat.ssl$U099_Count+1)+0.3055942*dat.ssl$model.score)


