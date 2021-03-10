# Generating survival times using Cox proportional hazards models with cyclic and piecewise time-varying covariates
This is the R code repository for implementing methods and simulations described in Huang Y, Zhang L, Zhang Z, Gilbert P. (2020) "Generating Survival Times Using Cox Proportional Hazards Models with Cyclic and Piecewise Time-Varying Covariates." Stat Biosci. 2020 Jan 25:1-16. doi: 10.1007/s12561-020-09266-3. Epub ahead of print. PMID: 32421033; PMCID: PMC7223425.

Specifically, it generates survival data with cyclic time-varying covariates using the antibody-mediated prevention (AMP) trials (ClinicalTrials.gov #NCT02716675 & #NCT02568215) as a concrete example, where HIV negative volunteers are randomized in 1:1:1 allocation to one of three study arms: 10 mg/kg VRC01 (low dose), 30 mg/kg VRC01 (high dose), or placebo, each administered via infusion every 8 weeks for 10 infusions. The cyclic time-varying covariate (e.g., time since last infusion in AMP) is first simulated before the survival time is simulated. 

For the simulations described in Section 3 of Huang et al. (2020), the following steps are carried out: 
1. Infusion visit times for each of the 10 infusions are simulated for each participant. More details about this step can be found in Zhang , Gilbert P, Capparelli E, Huang Y. arXiv:1801.08626, 2018.
2. The single-dose approach and the  multiple-dose approach described in Section 2.4 are used to simulate time to infection as a funtion of time since last infusion. 
3. The final survival data is taken as the minimal of the simulated time to infection and the end of follow up time (e.g, 80 weeks). 

## Authors: 
Lily Zhang, Yunda Huang
