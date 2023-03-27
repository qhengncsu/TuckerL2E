## Code from other sources
This repository also contains code from the following sources with our own modifications. The orignal sources and modifications are documented as following: 
- [HoRPCA](https://onedrive.live.com/?authkey=%21AOPu2g59n7NqZBI&id=731BCE806DD1BE58%216666&cid=731BCE806DD1BE58&parId=root&parQt=sharedby&o=OneUp) (Goldfarb, Donald, and Zhiwei Qin. "Robust low-rank tensor recovery: Models and algorithms." SIAM Journal on Matrix Analysis and Applications 35.1 (2014): 225-253.)

We have written several customized wrappers "twg_horpcas.m", "twg_rpca.m", "wrapper_horpcac.m", "wrapper_horpcas.m", "wrapper_rpca.m" to make it easier to invoke the respective methods.

- [BRTF](https://github.com/qbzhao/BRTF) (Zhao, Qibin, et al. "Bayesian robust tensor factorization for incomplete multiway data." IEEE transactions on neural networks and learning systems 27.4 (2015): 736-748.)

The videos in [BRTF](https://github.com/qbzhao/BRTF) are not included to save space. 

- [RGrad](https://www.tandfonline.com/doi/suppl/10.1080/01621459.2022.2063131?scroll=top&role=tab) (Cai, Jian-Feng, Jingyang Li, and Dong Xia. "Generalized low-rank plus sparse tensor estimation by fast Riemannian optimization." Journal of the American Statistical Association (2022): 1-17.)

We placed a function "wshift.m" in the "functions" subfolder to eliminate the dependency on the Wavelet Toolbox. We have renamed the "threshold" funciton to "mythreshold" to avoid a name conflict with the Econometrics Toolbox. The syntax for the size function in the function "tprod" is also changed to make the code excutable on earlier versions of Matlab.

- [tRPCA](https://github.com/canyilu/Tensor-Robust-Principal-Component-Analysis-TRPCA) (Lu, Canyi, et al. "Tensor robust principal component analysis with a new tensor nuclear norm." IEEE transactions on pattern analysis and machine intelligence 42.4 (2019): 925-938.)

"tprod.m" is renamed to "tprod_trpca.m" to avoid a name conflict with the software of [RGrad](https://www.tandfonline.com/doi/suppl/10.1080/01621459.2022.2063131?scroll=top&role=tab).
