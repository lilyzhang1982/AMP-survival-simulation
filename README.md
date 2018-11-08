# Generating survival times using Cox proportional hazards models with cyclic and piecewise time-varying covariates
This is the R code repository for implementing methods and simulations described in Huang Y, Zhang L, Zhang Z, Gilbert P. (2018) " Generating survival times using Cox proportional hazards models with cyclic time-varying covariates." arXiv:1801.08248, 2018.

Specifically, it generates survival data with cyclic time-varying covariates using the antibody-mediated prevention (AMP) trials as an example, where HIV negative volunteers are randomized in 1:1:1 allocation to one of three study arms{ 10 mg=kg VRC01 (low dose), 30 mg=kg VRC01 (high dose), or placebo, each administered via infusion every 8 weeks for 10 infusions. 

Consquently, the cyclic time-varying covariate (e.g., time since last infusion in AMP) is first simulated before the survival time is simulated. For the simulations described in Section 3 of Huang et al. (2018), The following steps are carried out: AMP example, time-since-last infusion 
1. Infusion visit times for each of the 10 infusions are simulated for each participant assuming perfect study product adherence
2. The single-dose approach described in Section 2.4.1 is used to simulate time to infection as a funtion of time since last infusion. 
3. The final survival data is taken as the minimal of the simulated time to infection and 80 weeks. 

## Author: 
Lily Zhang, Yunda Huang
