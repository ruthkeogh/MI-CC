

# ------------------------------
# packages

library(survival)
library(mice)
library(smcfcs)
library(mitools)

# ------------------------------
# number of imputations

nimp=10

#===============================
#complete case analysis using case-cohort data
#===============================

# load data set
load(file="caco.RData") #this is 'caco'

# fit the model
model=coxph(Surv(entertime,t,d)~x1+x2+z+cluster(id),data=caco)
summary(model)

#===============================
#MI-approx:  full-cohort approach
#===============================

# load data set
load(file="cohort_caco.RData") #this is 'cohort.caco'

# Compute Nelson-Aalen estimate of the cumulative hazard
cohort.caco$chaz=nelsonaalen(cohort.caco,t,d)

#predictor matrix which determines the imputation models for x1 and x2
pred.mat=matrix(0,nrow=dim(cohort.caco)[2],ncol=dim(cohort.caco)[2])
colnames(pred.mat)=names(cohort.caco)
rownames(pred.mat)=names(cohort.caco)
pred.mat["x1",c("x2","z","d","chaz")]=1
pred.mat["x2",c("x1","z","d","chaz")]=1

#method of imputation for x1 and x2
method.vec=rep("",dim(cohort.caco)[2])
method.vec[which(colnames(cohort.caco)=="x1")]="norm"
method.vec[which(colnames(cohort.caco)=="x2")]="logreg"

#make x2 a factor for the imputation
cohort.caco$x2<-as.factor(cohort.caco$x2)

#perform the imputation
imp<-mice(cohort.caco, m = nimp, method = method.vec, 
          predictorMatrix = pred.mat, 
          maxit = 5, diagnostics = FALSE, printFlag = T)

# Fit the analysis model in each imputed data set
models<-with(imp,coxph(Surv(t,d)~x1+x2+z))

# Combine estimates across the imputed data sets using Rubin's Rules
summary(pool(models,method="plain"))

#===============================
#MI-Approx: intermediate approach
#===============================

# obtain estimates across the imputed data sets and combine using Rubin's Rules
models<- with(imp, coxph(Surv(entertime[subco==1|d==1],t[in.caco==1],d[in.caco==1],type="counting")~
                               x1[in.caco==1]+x2[in.caco==1]+z[in.caco==1]+cluster(id[in.caco==1])))
summary(pool(models,method="plain"))

#===============================
#MI-SMC: full-cohort approach
#===============================

# load data set
load(file="cohort_caco.RData") #this is 'cohort.caco'

#predictor matrix which determines the imputation models for x1 and x2 (the outcomes are not included - see smcfcs help)
pred.mat=matrix(0,nrow=dim(cohort.caco)[2],ncol=dim(cohort.caco)[2])
colnames(pred.mat)=names(cohort.caco)
rownames(pred.mat)=names(cohort.caco)
pred.mat["x1",c("x2","z")]=1
pred.mat["x2",c("x1","z")]=1

#method of imputation for x1 and x2
method.vec=rep("",dim(cohort.caco)[2])
method.vec[which(colnames(cohort.caco)=="x1")]="norm"
method.vec[which(colnames(cohort.caco)=="x2")]="logreg"

#make x2 a factor for the imputation
cohort.caco$x2<-as.factor(cohort.caco$x2)

#perform the imputation
imp <- smcfcs(cohort.caco, smtype="coxph", smformula="Surv(t,d)~x1+x2+z",method=method.vec,
              predictorMatrix=pred.mat,m = nimp, numit = 100, rjlimit = 10000,noisy=F)

# obtain estimates across the imputed data sets and combine using Rubin's Rules
impobj <- imputationList(imp$impDatasets)
models <- with(impobj, coxph(Surv(t,d)~x1+x2+z))
MIcombine(models)

#===============================
#MI-SMC: intermediate approach
#===============================

# obtain estimates across the imputed data sets and combine using Rubin's Rules
impobj <- imputationList(imp$impDatasets)
models <- with(impobj, coxph(Surv(entertime[in.caco==1],t[in.caco==1],d[in.caco==1],type="counting")~
                               x1[in.caco==1]+x2[in.caco==1]+z[in.caco==1]+cluster(id[in.caco==1])))
MIcombine(models)

#===============================
#MI-approx: substudy approach
#===============================

# load data set
load(file="caco.RData") #this is 'caco'

# Compute Nelson-Aalen-type estimate H*{CC} (see paper section 4.3.1)
caco$chaz=nelsonaalen(caco,t,d)

#predictor matrix which determines the imputation model for x2
pred.mat=matrix(0,nrow=dim(caco)[2],ncol=dim(caco)[2])
colnames(pred.mat)=names(caco)
rownames(pred.mat)=names(caco)
pred.mat["x2",c("x1","z","d","chaz")]=1

#method of imputation for x2
method.vec=rep("",dim(caco)[2])
method.vec[which(colnames(caco)=="x2")]="logreg"

#make x2 a factor for the imputation
caco$x2<-as.factor(caco$x2)

#perform the imputation
imp<-mice(caco, m = nimp, method = method.vec, 
          predictorMatrix = pred.mat, 
          maxit = 5, diagnostics = FALSE, printFlag = T)

# Fit the analysis model in each imputed data set
models<-with(imp,coxph(Surv(entertime,t,d)~x1+x2+z+cluster(id)))

# Combine estimates across the imputed data sets using Rubin's Rules
summary(pool(models,method="plain"))

#===============================
#MI-SMC: substudy approach
#===============================

# load data set
load(file="caco.RData") #this is 'caco'

#proportion of the full cohort who are in the subcohort
my.sampfrac<-750/15000 

#predictor matrix which determines the imputation model for x2
pred.mat=matrix(0,nrow=dim(caco)[2],ncol=dim(caco)[2])
colnames(pred.mat)=names(caco)
rownames(pred.mat)=names(caco)
pred.mat["x2",c("x1","z")]=1

#method of imputation for x2
method.vec=rep("",dim(caco)[2])
method.vec[which(colnames(caco)=="x2")]="logreg"

#make x2 a factor for the imputation
caco$x2<-as.factor(caco$x2)

#perform the imputation
imp <- smcfcs.casecohort(caco,smformula="Surv(entertime,t,d)~x1+x2+z",sampfrac=my.sampfrac,in.subco="subco",
                         method=method.vec,predictorMatrix=pred.mat,m=nimp,numit=100,rjlimit=10000,noisy=FALSE)

# obtain estimates across the imputed data sets and combine using Rubin's Rules
impobj <- imputationList(imp$impDatasets)
models <- with(impobj, coxph(Surv(entertime,t,d)~x1+x2+z+cluster(id)))
MIcombine(models)
