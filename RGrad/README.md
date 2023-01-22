## Codes for  "Generalized Low-rank plus Sparse Tensor Estimation by Fast Riemannian Optimization" by Cai, Li and Xia, 2021

These codes are for the synthetic data presented in Section 6 and Section 11.

The file `compare_pca.m` is for the first experiment in Section 6, where we compare RGrad and PGD in terms of the tensor PCA problem.

The file `compare_RGrad_PGD.m`, `compare_hooi.m` are for the comparison of our method with other methods over the tensor RPCA problem.

The file `compare_binary.m` is for the comparison of our method with PGD over the binary tensor recovery problem.

The file `test_tensorRPCA.m` is for the experiments presented in Section 11.2 where we consider the robust Tensor PCA.

The file `test_heavy_tail_tensorRPCA.m` is for the experiments presented in Section 11.3 where we consider the Tensor PCA with heavy tail noise.

The file `test_binary.m` is for the experiments presented in Section 11.4 where we consider the binary tensor.

The file `test_poisson_rpca.m` is for the experiments presented in Section 11.5 where we consider the Poisson tensor RPCA.





These codes are for the BIC test on synthetic data and real data.

The files `BICtest_rpca.m` and `BICtest_binary.m` are for the Section 11.1.

The files `BICtest_realdata.m` and `test_real_data.m` are for the Section 12.

The files `BIC_compare_rpca.m` and `BIC_compare_binary.m` are for the Section 6, where we use BIC to select the parameters ($\alpha$ and $r$).





These codes are for the real data analysis.

The file `tensorRPCA.m` in the folder Real Data Analysis/International Trade is used for analyzing the internaction trade flow dataset and draw the plots

The file `BICtest_realdata.m` in the folder Real Data Analysis/International Trade is used for evaluating the BIC values for the international trade flow dataset

The file `test_real_data.m` in the folder Real Data Analysis/International Trade is used for evaluating the prediction performance for the international trade flow dataset

The file `tensorRPCA_SCORE.m` in the folder Real Data Analysis/Statistician Network is used for analyzing the statistician coauthorship network and draw the plots


!!!!!!!!!!!!!!   Remember to addpath("functions") before running the above files
