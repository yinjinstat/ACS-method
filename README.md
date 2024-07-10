# ACS-method
Code for Paper "A unified generalization of inverse regression via adaptive column selection"

newest_ACS_SDR.R includes the updated code used for the main experiments, the main function is named as "Four_estimate_c"(But only includes the ACS-method applied to M_{SAVE} and (M_{SAVE},M_{SIR})

SEAS.r: Jing Zeng's code for SEAS approach in paper "Subspace estimation with automatic dimension and variable selection in sufficient dimension reduction", can be also found in his github page https://github.com/jingzzeng/SEAS.

GSE24417_Training_Datamatrix.csv: the predictors of the dataset used in the real data analysis.

trainy.csv: the response of the dataset used in the real data analysis.

SC-SIR_estimate.csv: the estimate of central subspace produced by SC-SIR method from its matlab code.

perfect_outcome.csv: the estimate of central subspace produced by our AS method, which was also used to generate the graph in the paper.

SEAS-SIR_estimate.csv the estimate of central subspace produced by SEAS-SIR method 

tcsavegamma.csv the estiamte of central subspace produced by TC-sAVE method from its matlab code.
