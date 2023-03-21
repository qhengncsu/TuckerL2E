# Code for  "Robust Low-rank Tensor Decomposition with the L2 Criterion"

Qiang Heng, Eric C. Chi, Yufeng Liu (2023).
This repository is maintained by Qiang Heng (qheng@ncsu.edu).

## Dependencies
- [Matlab 2019b](https://www.mathworks.com/products/matlab.html) or newer
- [Statistics and Machine Learning Toolbox](https://www.mathworks.com/products/statistics.html), which needs to be installed with the Add-on Explorer within Matlab.
- This repo is also dependent on [Tensor Toolbox for MATLAB v3.2.1](https://www.tensortoolbox.org/) and the [L-BFGS-B implementation by Stephen Becker](https://github.com/stephenbeckr/L-BFGS-B-C), which has been included in "[folder name]/TuckerL2E" (no action is needed on your end). The "hosvd" function in the tensor toolbox has been renamed to "hosvd_tt" to avoid a name conflict.

## Compiling the mex files
After downloading and unzipping this repo, you will also need to compile the mex files for [L-BFGS-B](https://github.com/stephenbeckr/L-BFGS-B-C) before you attempt to run the scripts.

1. For Windows users, please install "MATLAB Support for MinGW-w64 C/C++ Compiler" using the Add-on Explorer within Matlab.
2. For Mac users, please ensure Xcode is installed and functional.
3. For Linux users, our experience has been that the binary files provided are typically directly usable and you can skip compilation. But if that is not the case, make sure GCC compiler is installed (most linux distributions have GCC pre-installed) and go to step 4.
4. Navigate to "[folder name]/TuckerL2E/L-BFGS-B-C/Matlab" and run the script "compile_mex.m" within Matlab (The working directory of Matlab should be "[folder name]/TuckerL2E/L-BFGS-B-C/Matlab"). 

## Example Scripts
Now you are ready to try out the scripts!

1. "generate_figure_2.m" will help you recreate the fourth subplot of Figure 2 in the paper (using only one random replicate). 
2. "generate_figure_6.m" will help you recreate Figure 6 in the paper. 
3. "generate_figure_8.m" will help you recreate Figure 8 in the paper. 

"generate_figure_8.m" will run relatively fast, "generate_figure_6.m" will take about 10-15 minutes, "generate_figure_2.m" will take about 15-25 minutes.

## Code from other sources
This repository also contains code from the following sources with some of our own modifications. The orignal sources and modifications are documented as following: 
- [HoRPCA](https://onedrive.live.com/?authkey=%21AOPu2g59n7NqZBI&id=731BCE806DD1BE58%216666&cid=731BCE806DD1BE58&parId=root&parQt=sharedby&o=OneUp) (Goldfarb, Donald, and Zhiwei Qin. "Robust low-rank tensor recovery: Models and algorithms." SIAM Journal on Matrix Analysis and Applications 35.1 (2014): 225-253.)

We have written several customized wrappers "twg_horpcas.m", "twg_rpca.m", "wrapper_horpcac.m", "wrapper_horpcas.m", "wrapper_rpca.m" to make it easier to invoke the respective methods.

- [RGrad](https://www.tandfonline.com/doi/suppl/10.1080/01621459.2022.2063131?scroll=top&role=tab) (Cai, Jian-Feng, Jingyang Li, and Dong Xia. "Generalized low-rank plus sparse tensor estimation by fast Riemannian optimization." Journal of the American Statistical Association (2022): 1-17.)

We placed a function "wshift.m" in the "functions" subfolder to eliminate the dependency on the Signal Processing Toolbox. We have renamed the "threshold" funciton to "mythreshold" to avoid a name conflict with the Econometrics Toolbox. The syntax for the size function in the function "tprod" is also changed to make the code excutable on earlier versions of Matlab.

- [tRPCA](https://github.com/canyilu/Tensor-Robust-Principal-Component-Analysis-TRPCA) (Lu, Canyi, et al. "Tensor robust principal component analysis with a new tensor nuclear norm." IEEE transactions on pattern analysis and machine intelligence 42.4 (2019): 925-938.)
