

# ------------------------------
# packages

library(survival)
library(mice)
library(smcfcs)
library(mitools)
library(multipleNCC)

# ------------------------------
# number of imputations

nimp=10

#===============================
#complete case analysis using nested case-control data
#===============================

# load data set
load(file="ncc.RData") #this is 'ncc'

# fit the model
model=coxph(Surv(t,case)~x1+x2+z+strata(setno),data=ncc)
summary(model)

#===============================
#MI-approx:  full-cohort approach
#===============================

# load data set
load(file="cohort_ncc.RData") #this is 'cohort.ncc'

# Compute Nelson-Aalen estimate of the cumulative hazard
cohort.ncc$chaz=nelsonaalen(cohort.ncc,t,d)

#predictor matrix which determines the imputation models for x1 and x2
pred.mat=matrix(0,nrow=dim(cohort.ncc)[2],ncol=dim(cohort.ncc)[2])
colnames(pred.mat)=names(cohort.ncc)
rownames(pred.mat)=names(cohort.ncc)
pred.mat["x1",c("x2","z","d","chaz")]=1
pred.mat["x2",c("x1","z","d","chaz")]=1

#method of imputation for x1 and x2
method.vec=rep("",dim(cohort.ncc)[2])
method.vec[which(colnames(cohort.ncc)=="x1")]="norm"
method.vec[which(colnames(cohort.ncc)=="x2")]="logreg"

#make x2 a factor for the imputation
cohort.ncc$x2<-as.factor(cohort.ncc$x2)

#perform the imputation 
imp<-mice(cohort.ncc, m = nimp, method = method.vec, 
          predictorMatrix = pred.mat, 
          maxit = 5, diagnostics = FALSE, printFlag = T)

# Fit the analysis model in each imputed data set
model<-with(imp,coxph(Surv(t,d)~x1+x2+z))

# Combine estimates across the imputed data sets using Rubin's Rules
summary(pool(model,method="plain"))

#===============================
#MI-Approx: intermediate approach
#===============================

temp.ncc<-ncc[,c("id","t","d","x1","z","setno","case")]#this omits the variable which has been imputed

coef.imp.est<-matrix(nrow=nimp,ncol=3)
coef.imp.var<-matrix(nrow=nimp,ncol=3)

for (k in 1:nimp){
  temp.imp<-complete(imp,k)[,c("id","x2")]
  
  ncc.imp<-merge(temp.ncc,temp.imp,by.x="id")
  
  model=coxph(Surv(t,case)~x1+x2+z+strata(setno),data=ncc.imp)
  coef.imp.est[k,]<-model$coef
  coef.imp.var[k,]<-diag(model$var)
}

pool.scalar(coef.imp.est[,1],coef.imp.var[,1])
pool.scalar(coef.imp.est[,2],coef.imp.var[,2])
pool.scalar(coef.imp.est[,3],coef.imp.var[,3])

#===============================
#MI-SMC: full-cohort approach
#===============================

# load data set
load(file="cohort_ncc.RData") #this is 'cohort.ncc'

#predictor matrix which determines the imputation models for x1 and x2 (the outcomes are not included - see smcfcs help)
pred.mat=matrix(0,nrow=dim(cohort.ncc)[2],ncol=dim(cohort.ncc)[2])
colnames(pred.mat)=names(cohort.ncc)
rownames(pred.mat)=names(cohort.ncc)
pred.mat["x1",c("x2","z")]=1
pred.mat["x2",c("x1","z")]=1

#method of imputation for x1 and x2
method.vec=rep("",dim(cohort.ncc)[2])
method.vec[which(colnames(cohort.ncc)=="x1")]="norm"
method.vec[which(colnames(cohort.ncc)=="x2")]="logreg"

#make x2 a factor for the imputation
cohort.ncc$x2<-as.factor(cohort.ncc$x2)

#perform the imputation
imp <- smcfcs(cohort.ncc, smtype="coxph", smformula="Surv(t,d)~x1+x2+z",method=method.vec,
              predictorMatrix=pred.mat,m = nimp, numit = 100, rjlimit = 10000,noisy=F)

# obtain estimates across the imputed data sets and combine using Rubin's Rules
impobj <- imputationList(imp$impDatasets)
models <- with(impobj, coxph(Surv(t,d)~x1+x2+z))
MIcombine(models)

#===============================
#MI-SMC: intermediate approach
#===============================

temp.ncc<-ncc[,c("id","t","d","x1","z","setno","case")]#this omits x2

coef.imp.est<-matrix(nrow=nimp,ncol=3)
coef.imp.var<-matrix(nrow=nimp,ncol=3)

for (k in 1:nimp){
  temp.imp<-imp$impDatasets[[k]][,c("id","x2")]
  
  ncc.imp<-merge(temp.ncc,temp.imp,by.x="id")
  
  model=coxph(Surv(t,case)~x1+x2+z+strata(setno),data=ncc.imp)
  coef.imp.est[k,]<-model$coef
  coef.imp.var[k,]<-diag(model$var)
}

pool.scalar(coef.imp.est[,1],coef.imp.var[,1])
pool.scalar(coef.imp.est[,2],coef.imp.var[,2])
pool.scalar(coef.imp.est[,3],coef.imp.var[,3])

#===============================
#MI-approx: substudy approach
#this method uses the full cohort information on (T,D) to estimate the cumulative hazard
#===============================

# load data sets
load(file="cohort_ncc.RData") #this is 'cohort.ncc'
load(file="ncc.RData") #this is 'ncc'

# Compute Nelson-Aalen estimate of the cumulative baseline hazard
cohort.ncc$chaz=nelsonaalen(cohort.ncc,t,d)

#add cumulative hazard into ncc data
cohort.ncc.merge<-cohort.ncc[,c("id","chaz")] 
ncc<-merge(ncc,cohort.ncc.merge,by.x="id")

#predictor matrix which determines the imputation model for x2
pred.mat=matrix(0,nrow=dim(ncc)[2],ncol=dim(ncc)[2])
colnames(pred.mat)=names(ncc)
rownames(pred.mat)=names(ncc)
pred.mat["x1",c("x2","z","d","chaz")]=1
pred.mat["x2",c("x1","z","d","chaz")]=1

#method of imputation for x2
method.vec=rep("",dim(ncc)[2])
method.vec[which(colnames(ncc)=="x2")]="logreg"

#make x2 a factor for the imputation
ncc$x2<-as.factor(ncc$x2)

#perform the imputation
imp<-mice(ncc, m = nimp, method = method.vec, 
          predictorMatrix = pred.mat, 
          maxit = 5, diagnostics = FALSE, printFlag = T)

# Fit the analysis model in each imputed data set
model<-with(imp,coxph(Surv(t,case)~x1+x2+z+strata(setno)))

# Combine estimates across the imputed data sets using Rubin's Rules
summary(pool(model,method="plain"))

#===============================
#MI-SMC: substudy approach
#The model for the proposal distribution is fitted in the 'true' controls only
#===============================

# load data sets
load(file="cohort_ncc.RData") #this is 'cohort.ncc'
load(file="ncc.RData") #this is 'ncc'

# Compute number at risk at each event time using the full cohort data
nrisk.fit<-survfit(Surv(t,d)~1,data=cohort.ncc)
ord.t.d1<-order(cohort.ncc$t[cohort.ncc$d==1])
numrisk<-summary(nrisk.fit,censored=F)$n.risk #this gives the number at risk at each unique event time

#add numbers at risk at each event time into the nested case-control data set
ncc$numrisk<-NA
ncc$numrisk[ncc$case==1][ord.t.d1]<-numrisk

#In each matched set: assign number at risk at the case's event time to every individual in the set
ncc$numrisk<-ave(ncc$numrisk, ncc$setno, FUN = function(x) sum(x, na.rm=T))

#predictor matrix which determines the imputation model for x2
pred.mat=matrix(0,nrow=dim(ncc)[2],ncol=dim(ncc)[2])
colnames(pred.mat)=names(ncc)
rownames(pred.mat)=names(ncc)
pred.mat["x2",c("x1","z")]=1

#method of imputation for x2
method<-rep("",dim(ncc)[2])
method[which(colnames(ncc)=="x2")]="logreg"

#perform the imputation
imp<-smcfcs.nestedcc(ncc,smformula="Surv(t,case)~x1+x2+z+strata(setno)",set="setno",event="d",nrisk="numrisk",
                     method=method,predictorMatrix=pred.mat,m=nimp,numit=10,rjlimit=1000,noisy=F) 

impobj <- imputationList(imp$impDatasets)
models <- with(impobj, coxph(Surv(t,case)~x1+x2+z+strata(setno)))
MIcombine(models)

#=====================================================================================================
#MI using 'matched set'
#=====================================================================================================

# load data set
load(file="ncc.RData") #this is 'ncc'

#order the ncc data by setno and case

ncc<-ncc[order(ncc$setno,1-ncc$case),]

#drop id

ncc<-ncc[,-c(1)]

#generate a within-set id number which is 1 for the case in each set and 2-5 for the 4 controls

vec.ones<-rep(1,dim(ncc)[1])
ncc$setid<-ave(vec.ones, ncc$setno, FUN = cumsum)

#reshape the data in 'wide' form

ncc.wide<-reshape(ncc,v.names=c("t","d","x1","x2","z","case"),timevar="setid",idvar="setno",direction="wide")

#generate partial sums of x1, x2 and z in each set - each time excluding one member of the set
ncc.wide$x1sum.1<-ncc.wide$x1.2+ncc.wide$x1.3+ncc.wide$x1.4+ncc.wide$x1.5
ncc.wide$x1sum.2<-ncc.wide$x1.1+ncc.wide$x1.3+ncc.wide$x1.4+ncc.wide$x1.5
ncc.wide$x1sum.3<-ncc.wide$x1.1+ncc.wide$x1.2+ncc.wide$x1.4+ncc.wide$x1.5
ncc.wide$x1sum.4<-ncc.wide$x1.1+ncc.wide$x1.2+ncc.wide$x1.3+ncc.wide$x1.5
ncc.wide$x1sum.5<-ncc.wide$x1.1+ncc.wide$x1.2+ncc.wide$x1.3+ncc.wide$x1.4

ncc.wide$x2sum.1<-ncc.wide$x2.2+ncc.wide$x2.3+ncc.wide$x2.4+ncc.wide$x2.5
ncc.wide$x2sum.2<-ncc.wide$x2.1+ncc.wide$x2.3+ncc.wide$x2.4+ncc.wide$x2.5
ncc.wide$x2sum.3<-ncc.wide$x2.1+ncc.wide$x2.2+ncc.wide$x2.4+ncc.wide$x2.5
ncc.wide$x2sum.4<-ncc.wide$x2.1+ncc.wide$x2.2+ncc.wide$x2.3+ncc.wide$x2.5
ncc.wide$x2sum.5<-ncc.wide$x2.1+ncc.wide$x2.2+ncc.wide$x2.3+ncc.wide$x2.4

ncc.wide$zsum.1<-ncc.wide$z.2+ncc.wide$z.3+ncc.wide$z.4+ncc.wide$z.5
ncc.wide$zsum.2<-ncc.wide$z.1+ncc.wide$z.3+ncc.wide$z.4+ncc.wide$z.5
ncc.wide$zsum.3<-ncc.wide$z.1+ncc.wide$z.2+ncc.wide$z.4+ncc.wide$z.5
ncc.wide$zsum.4<-ncc.wide$z.1+ncc.wide$z.2+ncc.wide$z.3+ncc.wide$z.5
ncc.wide$zsum.5<-ncc.wide$z.1+ncc.wide$z.2+ncc.wide$z.3+ncc.wide$z.4

#predictor matrix which determines the imputation model for x2.k (k=1,2,3,4,5)
pred.mat<-matrix(0,nrow=dim(ncc.wide)[2],ncol=dim(ncc.wide)[2])
colnames(pred.mat)=names(ncc.wide)
rownames(pred.mat)=names(ncc.wide)
pred.mat["x2.1",c("x1.1","z.1","x1sum.1","x2sum.1","zsum.1")]<-1
pred.mat["x2.2",c("x1.2","z.2","x1sum.2","x2sum.2","zsum.2")]<-1
pred.mat["x2.3",c("x1.3","z.3","x1sum.3","x2sum.3","zsum.3")]<-1
pred.mat["x2.4",c("x1.4","z.4","x1sum.4","x2sum.4","zsum.4")]<-1
pred.mat["x2.5",c("x1.5","z.5","x1sum.5","x2sum.5","zsum.5")]<-1

#make x2.k (k=1,2,3,4,5) into factors for the imputation
ncc.wide$x2.1<-as.factor(ncc.wide$x2.1)
ncc.wide$x2.2<-as.factor(ncc.wide$x2.2)
ncc.wide$x2.3<-as.factor(ncc.wide$x2.3)
ncc.wide$x2.4<-as.factor(ncc.wide$x2.4)
ncc.wide$x2.5<-as.factor(ncc.wide$x2.5)

#method of imputation for x2.k (k=1,2,3,4,5) 
method.vec=rep("",dim(ncc.wide)[2])
method.vec[c(which(colnames(ncc.wide)=="x2.1"),which(colnames(ncc.wide)=="x2.2"),
             which(colnames(ncc.wide)=="x2.3"),which(colnames(ncc.wide)=="x2.4"),which(colnames(ncc.wide)=="x2.5"))]="logreg"

#set up the passive imputation of x2sum.k (k=1,2,3,4,5)
#note that the as.numeric elements as requried because the x2.k are factors
method.vec[which(colnames(ncc.wide)=="x2sum.1")]="~I(as.numeric(x2.2)-1+as.numeric(x2.3)-1+as.numeric(x2.4)-1+as.numeric(x2.5)-1)" 
method.vec[which(colnames(ncc.wide)=="x2sum.2")]="~I(as.numeric(x2.1)-1+as.numeric(x2.3)-1+as.numeric(x2.4)-1+as.numeric(x2.5)-1)" 
method.vec[which(colnames(ncc.wide)=="x2sum.3")]="~I(as.numeric(x2.1)-1+as.numeric(x2.2)-1+as.numeric(x2.4)-1+as.numeric(x2.5)-1)" 
method.vec[which(colnames(ncc.wide)=="x2sum.4")]="~I(as.numeric(x2.1)-1+as.numeric(x2.2)-1+as.numeric(x2.3)-1+as.numeric(x2.5)-1)" 
method.vec[which(colnames(ncc.wide)=="x2sum.5")]="~I(as.numeric(x2.1)-1+as.numeric(x2.2)-1+as.numeric(x2.3)-1+as.numeric(x2.4)-1)" 

imp<-mice(ncc.wide, m = nimp, method = method.vec, 
          predictorMatrix = pred.mat, 
          maxit = 5, diagnostics = FALSE, printFlag = T)

#reshape the imputed data sets into long format and fit the substantive model in each imputed data set

coef.imp.est<-matrix(nrow=nimp,ncol=3)
coef.imp.var<-matrix(nrow=nimp,ncol=3)

for (k in 1:nimp){
  temp.imp<-complete(imp,k)
  
  temp.imp.long<-reshape(temp.imp,varying= c("t.1","d.1","x1.1","z.1","x2.1","case.1","x1sum.1","x2sum.1","zsum.1",
                                             "t.2","d.2","x1.2","z.2","x2.2","case.2","x1sum.2","x2sum.2","zsum.2",
                                             "t.3","d.3","x1.3","z.3","x2.3","case.3","x1sum.3","x2sum.3","zsum.3",
                                             "t.4","d.4","x1.4","z.4","x2.4","case.4","x1sum.4","x2sum.4","zsum.4",
                                             "t.5","d.5","x1.5","z.5","x2.5","case.5","x1sum.5","x2sum.5","zsum.5"),
                         direction="long",idvar="setno",sep=".",timevar="setid")
  
  model=coxph(Surv(t,case)~x1+x2+z+strata(setno),data=temp.imp.long)
  coef.imp.est[k,]<-model$coef
  coef.imp.var[k,]<-diag(model$var)
}

# Combine estimates across the imputed data sets using Rubin's Rules
pool.scalar(coef.imp.est[,1],coef.imp.var[,1])
pool.scalar(coef.imp.est[,2],coef.imp.var[,2])
pool.scalar(coef.imp.est[,3],coef.imp.var[,3])

