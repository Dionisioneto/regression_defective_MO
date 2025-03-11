**Estimating cure rates with Marshall-Olkin defective regression models**

This repository contains code for the manuscript Estimating cure rates with Marshall-Olkin defective regression models by Dionisio Silva Neto, Vera Tomazella and Adriano Suzuki, in the process of submission.



📁**MO_Gompertz**

MLE_reg_MO_Gompertz.R contains the code to conduct the maximum likelihood estimation (MLE) of the Marshall-Olkin Gompertz regression model.

Monte_Carlo_frequentist_reg_MO_Gompertz.R contains the iterative code to conduct the Monte carlos simulation of the Marshall-Olkin Gompertz regression model.

resumos_freq_MC_reg_MO_Gompertz is a folder which presents the results for different sample sizes of 1000 Monte Carlo replicates in the simulation study.


📁**MO_IG**

MLE_reg_MO_IG.R contains the code to conduct the maximum likelihood estimation (MLE) of the Marshall-Olkin Inverse Gaussian regression model.

Monte_Carlo_frequentist_reg_MO_IG.R contains the iterative code to conduct the Monte carlos simulation of the Marshall-Olkin Inverse Gaussian regression model.

resumos_freq_MC_reg_MO_IG is a folder which presents the results for different sample sizes of 1000 Monte Carlo replicates in the simulation study.


📁**application**

application1_freq_Colon_cancer.R contains the code necessary to replicate the data analysis in colon cancer data, using the defective regression model developed in the paper, the classical defective regression models
and the mixture cure weibull model. 

compare_models_Colon2.xlsx contains the table of consistent AIC, BIC, CAIC and others comparative measures of each model fitted to colon data.

reestimation_global_influence.xlsx  contains the summary of relative changes in the estimation of the chosen model for colon data (Marshall-Olkin Gompertz regression) in global influence studies.

source_reg_MO_Gompertz.R contains the source functions for estimation of Marshall-Olkin Gompertz regression.

source_reg_MO_IG.R contains the source functions for estimation of Marshall-Olkin inverse Gaussian regression.

source_reg_classical_defective.R contains the source functions for estimation of Gompertz and inverse Gaussian regression models.

Note: The dataset considered for this study is fully available in the survival package.
