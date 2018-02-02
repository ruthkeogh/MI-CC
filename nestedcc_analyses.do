
*=====================================================================================================
*complete case analysis using nested case-control data
*=====================================================================================================

* load data set and save as Stata data set for later methods
import delimited ncc.csv, numericcol(5) clear

stset t, failure(case)

stcox x1 x2 z, strata(setno) nohr

*=====================================================================================================
*MI-approx full cohort: imputation and analysis performed using full cohort
*=====================================================================================================

* load data set
import delimited cohort_ncc.csv, numericcol(4 5) clear

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
import delimited cohort_ncc.csv, numericcol(4 5) clear

stset t, failure(d)

*perform the imputation (10 imputations)
* Fit the analysis model in each imputed data set and combine estimates across the imputed data sets using Rubin's Rules
smcfcs stcox x1 x2 z, regress(x1) logit(x2) m(10) iterations(100)

*=====================================================================================================
*MI-approx nested case-control: imputation and analysis performed using nested case-control data only
*this method uses the full cohort to estimate the cumulative hazard
*x1 is fully observed in the NCC data
*=====================================================================================================

* load ncc data set and save - for use in merging below
import delimited ncc.csv, numericcol(4 5) clear

save ncc, replace

* load cohort data set
import delimited cohort_ncc.csv, numericcol(4 5) clear

stset t, failure(d)

* Compute Nelson-Aalen estimate of the cumulative baseline hazard
sts generate H=na

*add cumulative hazard into ncc data
merge 1:m id using ncc.dta

keep if _merge==3
drop _merge

*stset the data
stset t, failure(case)

*perform the imputation (10 imputations)

mi set mlong
mi register imputed x1 x2
mi register regular d t H z

mi impute chained (regress) x1 (logit) x2 = d H z, add(10)

* Fit the analysis model in each imputed data set and combine estimates across the imputed data sets using Rubin's Rules
mim: stcox x1 x2 z, strata(setno) nohr

*=====================================================================================================
*MI using 'matched set'
*imputation and analysis performed using nested case-control data only
*x1 is fully observed in the NCC data
*=====================================================================================================

* load data set
import delimited ncc.csv, numericcol(4 5) clear

*order the ncc data by setno and case
gen mcase=1-case
sort setno mcase
drop mcase

*drop id

drop id

*generate a within-set id number which is 1 for the case in each set and 2-5 for the 4 controls

by setno: gen setid=_n

*reshape the data in 'wide' form

reshape wide x1 x2 z t d case, i(setno) j(setid)

*generate partial sums of x1, x2 and z in each set - each time excluding one member of the set

foreach var in x1 x2 z{
gen `var'sum1=`var'2+`var'3+`var'4+`var'5
gen `var'sum2=`var'1+`var'3+`var'4+`var'5
gen `var'sum3=`var'1+`var'2+`var'4+`var'5
gen `var'sum4=`var'1+`var'2+`var'3+`var'5
gen `var'sum5=`var'1+`var'2+`var'3+`var'4
}

*perform the imputation

ice d1 x11 x21 z1 x1sum1 x2sum1 zsum1 /*
*/  d2 x12 x22 z2 x1sum2 x2sum2 zsum2 /*
*/  d3 x13 x23 z3 x1sum3 x2sum3 zsum3 /*
*/  d4 x14 x24 z4 x1sum4 x2sum4 zsum4 /*
*/  d5 x15 x25 z5 x1sum5 x2sum5 zsum5, /*
*/ eq(x21: x11 z1 x1sum1 x2sum1 zsum1,/* 
*/ x22: x12 z2 x1sum2 x2sum2 zsum2,/* 
*/ x23: x13 z3 x1sum3 x2sum3 zsum3,/* 
*/ x24: x14 z4 x1sum4 x2sum4 zsum4,/* 
*/ x25: x15 z5 x1sum5 x2sum5 zsum5) /*
*/ passive(x2sum1:(x22+x23+x24+x25)\/*
*/		   x2sum2:(x21+x23+x24+x25)\/*
*/		   x2sum3:(x21+x22+x24+x25)\/*
*/		   x2sum4:(x21+x22+x23+x25)\/*
*/		   x2sum5:(x21+x22+x23+x24)) m(10) saving(nccmi.dta, replace)

*reshape the imputed data sets into long format

use nccmi.dta, clear

quietly: mi import ice, automatic

quietly: mi reshape long x1 x2 z d t case, i(setno) j(setid)

* Fit the analysis model in each imputed data set and combine estimates across the imputed data sets using Rubin's Rules
mi stset t,failure(case)

mim: stcox x1 x2 z, strata(setno) nohr
