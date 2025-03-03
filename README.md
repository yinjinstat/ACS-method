# ACS-method
Code for Paper "A unified generalization of inverse regression methods via column selection"

________Main code___________________________________________________________________________________
newest_ACS_SDR.R includes the updated code used for the main experiments, the main function is named as "Four_estimate_c"(But only includes the ACS-method applied to M_{SAVE} and (M_{SAVE},M_{SIR})
SEAS.r: Jing Zeng's code for SEAS approach in paper "Subspace estimation with automatic dimension and variable selection in sufficient dimension reduction", can be also found in his github page https://github.com/jingzzeng/SEAS. 

--------The following files include the data, supplementary code, and the results used in the real data analysis in the paper-------
trainx.zip the predictors of the dataset used in the real data analysis. The transformed data from GSE24417_Training_Datamatrix, which could be downloaded at https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE24417. And see the transformation procedure in the manuscript
trainy.csv: the response of the dataset used in the real data analysis.
SCS_SDR_outcome.csv: the estimate of central subspace produced by our AS method, which was also used to generate the graph in the paper.
SEAS-SIR_estimate.csv the estimate of central subspace produced by SEAS-SIR method (Zeng jing et al. 2024)
SC-SIR_estimate.csv: the estimate of central subspace produced by SC-SIR (Wei qian et al. 2019) method 
tcsavegamma.csv the estiamte of central subspace produced by TC-sAVE method 
lasso_sir_outcome.csv the estiamte of central subspace produced by LASSO-SIR method (Qian lin et al. 2019) 
