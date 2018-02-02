*=====================================================================================================
*complete case analysis using case-cohort data
*=====================================================================================================

* load data set
import delimited caco.csv, numericcol(5) clear

stset t, failure(d) enter(entertime)

* fit the model

stcox x1 x2 z, robust nohr

*=====================================================================================================
*MI-approx full cohort: imputation and analysis performed using full cohort
*=====================================================================================================

* load data set
import delimited cohort_caco.csv, numericcol(4 5) clear

stset t, failure(d)

* Compute Nelson-Aalen estimate of the cumulative baseline hazard
sts generate H=na

*perform the imputation (10 imputations)

mi set mlong
mi register imputed x1 x2
mi register regular d t H z

mi impute chained (regress) x1 (logit) x2 = d H z, add(10)

* Fit the analysis model in each imputed data set and combine estimates across the imputed data sets using Rubin's Rules
mim: stcox x1 x2 z, nohr

*=====================================================================================================
*MI-SMC full cohort: imputation and analysis performed using full cohort
*=====================================================================================================

* load data set
import delimited cohort_caco.csv, numericcol(4 5) clear

stset t, failure(d)

*perform the imputation (10 imputations)
* Fit the analysis model in each imputed data set and combine estimates across the imputed data sets using Rubin's Rules
smcfcs stcox x1 x2 z, regress(x1) logit(x2) m(10) iterations(100)

*=====================================================================================================
*MI-approx case-cohort: imputation and analysis performed using case-cohort data 
*x1 is fully observed in the case-cohort data
*=====================================================================================================

* load data set
import delimited caco.csv, numericcol(5) clear

stset t, failure(d) enter(entertime)

* Compute Nelson-Aalen estimate of the cumulative baseline hazard
sts generate H=na

*perform the imputation (10 imputations)

mi set mlong
mi register imputed x1 x2
mi register regular d t H z

mi impute chained (regress) x1 (logit) x2 = d H z, add(10)

* Fit the analysis model in each imputed data set and combine estimates across the imputed data sets using Rubin's Rules
mim: stcox x1 x2 z, robust nohr

