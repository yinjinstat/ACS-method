# ACS-method
Code for Paper "A unified generalization of inverse regression via adaptive column selection"

ACS_prepare.r: some preliminary functions such as calculating the power of a positivedefinite matrix are built.

ctsave.r: codes for testing AE approach in low dimensional setting.

kernel_matrix_prepare.r: since we transform the penalized least square into a pseudo-lasso-problem, we prepare the pesudo-response for each common inverse regression approach.

ACS_final.r: the main functions of conducting AS-method in HD setting with order determination.

Experiments_simulated.r: Codes for building models and recording the result of AS-method and other comparisons in HD setting.

SEAS.r: Jing Zeng's code for SEAS approach in paper "Subspace estimation with automatic dimension and variable selection in sufficient dimension reduction", can be also found in his github page https://github.com/jingzzeng/SEAS.

GSE24417_Training_Datamatrix.csv: the predictors of the dataset used in the real data analysis.

trainy.csv: the response of the dataset used in the real data analysis.

SC-SIR_estimate.csv: the estimate of central subspace produced by SC-SIR method from its matlab code.

perfect_outcome.csv: the estimate of central subspace produced by our AS-Hybrid method, which was also used to generate the graph in the paper.


