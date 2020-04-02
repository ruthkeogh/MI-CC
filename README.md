# MI-CC
Multiple imputation (MI) for case-cohort and nested case-control studies

Example code in R and Stata are provided corresponding to the methods described in the paper:

Keogh RH, Seaman SR, Bartlett JW, Wood AM. Multiple imputation of missing data in nested case-control and case-cohort studies. Biometrics 2018; 74: 1438-1449. https://onlinelibrary.wiley.com/doi/full/10.1111/biom.12910, https://doi.org/10.1111/biom.12910

The file *generate_data* contains R code to simulate a single full cohort data set and case-cohort and nested case-control studies within it. The following data sets are created:

    cohort_caco: the cohort data with indicators of being in the subcohort and being in the case-cohort sample
    
    caco: the case-cohort substudy data 
    
    cohort_ncc: the cohort data with indicators of being in the nested case-control sample
    
    ncc: the nested case-control substudy data
    
The file *casecohort_analyses.R* illustrates MI analyses using the full-cohort approach, intermediate approach and substudy approach for case-cohort studies using both MI-Approx and MI-SMC.

The file *nestedcc_analyses.R* illustrates MI analyses using the full-cohort approach, intermediate approach and substudy approach for nested case-control studies using both MI-Approx and MI-SMC. Use of the MI matched set method is also illustrated for the substudy approach.

The corresponding files *casecohort_analyses.do* and *nestedcc_analyses.do* illustrate how to perform the same methods using Stata, with the exception of the substudy approach using MI-SMC.
