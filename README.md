# Code for  "Robust Low-rank Tensor Decomposition with the L2 Criterion"

Qiang Heng, Eric C. Chi, Yufeng Liu (2023).
This repository is maintained by Qiang Heng (qheng@ncsu.edu).

## Dependencies
- [Matlab 2019b](https://www.mathworks.com/products/matlab.html) or newer
- [Statistics and Machine Learning Toolbox](https://www.mathworks.com/products/statistics.html), which needs to be installed with the Add-on Explorer within Matlab.

Steps to get the code up and running:
1. Install [Tensor Toolbox for MATLAB v3.2.1](https://gitlab.com/tensors/tensor_toolbox/-/releases/v3.2.1). Please download the zip file [Source code (zip)](https://gitlab.com/tensors/tensor_toolbox/-/archive/v3.2.1/tensor_toolbox-v3.2.1.zip) and unzip it. Then please place the folder at "[path of repo]/TuckerL2E/tensor_toolbox-v3.2.1". An example absolute path on my machine is "C:/Users/Qiang/Documents/TuckerL2E/TuckerL2E/tensor_toolbox-v3.2.1", where "C:/Users/Qiang/Documents/TuckerL2E" is the path for this repo. 
2. Please rename the "hosvd.m" function inside "[path of repo]/TuckerL2E/tensor_toolbox-v3.2.1" to "hosvd_tt.m". The function name written inside "hosvd_tt.m" must be changed as well. In other words, the first line of "hosvd_tt.m" should now be "function T = hosvd_tt(X,tol,varargin)".
3. Download as zip the github repo [L-BFGS-B-C](https://github.com/stephenbeckr/L-BFGS-B-C). Unzip the zip file. Notice that the folder name will be "L-BFGS-B-C-master" after being downloaded from github. Please rename the folder to "L-BFGS-B-C" and place the folder at "[path of repo]/TuckerL2E/L-BFGS-B-C". An example absolute path on my machine is "C:/Users/Qiang/Documents/TuckerL2E/TuckerL2E/L-BFGS-B-C".

## Compiling the mex files
The [L-BFGS-B-C](https://github.com/stephenbeckr/L-BFGS-B-C) repo requires compilation! Thus you must follow the following steps to compile the mex files:

4. For Windows users, please install "MATLAB Support for MinGW-w64 C/C++ Compiler" using the Add-on Explorer within Matlab.
5. For Mac users, please ensure Xcode is installed and functional.
6. For Linux users, make sure GCC compiler is installed (most linux distributions have GCC pre-installed).
7. Navigate to "[path of repo]/TuckerL2E/L-BFGS-B-C/Matlab" and run the script "compile_mex.m" within Matlab (The working directory of Matlab should be "[path of repo]/TuckerL2E/L-BFGS-B-C/Matlab"). 

## Example Scripts
Now you are ready to try out the scripts! (Remember to swith the working directory back to [path of repo] after step 7)

8. "generate_figure_2_subplot_[i].m" (i=1,2,3,4) will help you recreate the four subplots of Figure 2 in the paper. 
9. "generate_figure_6.m" will help you recreate Figure 6 in the paper. 
10. "generate_figure_8.m" will help you recreate Figure 8 in the paper. 

"generate_figure_8.m" will run relatively fast, "generate_figure_6.m" will take about 10-15 minutes, each one of "generate_figure_2_subplot_[i].m" will take about 15-25 minutes if the number of replicate "N" is set to 1. Notice that in the paper, Figure 2 is generated with 50 random replicates. However, one random replicate is enough to generate a figure that is similar to Figure 2. 

## Code from other sources
This repository also contains code from the following sources with some of our own modifications. The orignal sources and modifications are documented as following: 
- [HoRPCA](https://onedrive.live.com/?authkey=%21AOPu2g59n7NqZBI&id=731BCE806DD1BE58%216666&cid=731BCE806DD1BE58&parId=root&parQt=sharedby&o=OneUp) (Goldfarb, Donald, and Zhiwei Qin. "Robust low-rank tensor recovery: Models and algorithms." SIAM Journal on Matrix Analysis and Applications 35.1 (2014): 225-253.)

We have written several customized wrappers "twg_horpcas.m", "twg_rpca.m", "wrapper_horpcac.m", "wrapper_horpcas.m", "wrapper_rpca.m" to make it easier to invoke the respective methods.

- [BRTF](https://github.com/qbzhao/BRTF) (Zhao, Qibin, et al. "Bayesian robust tensor factorization for incomplete multiway data." IEEE transactions on neural networks and learning systems 27.4 (2015): 736-748.)

- [RGrad](https://www.tandfonline.com/doi/suppl/10.1080/01621459.2022.2063131?scroll=top&role=tab) (Cai, Jian-Feng, Jingyang Li, and Dong Xia. "Generalized low-rank plus sparse tensor estimation by fast Riemannian optimization." Journal of the American Statistical Association (2022): 1-17.)

We placed a function "wshift.m" in the "functions" subfolder to eliminate the dependency on the Wavelet Toolbox. We have renamed the "threshold" funciton to "mythreshold" to avoid a name conflict with the Econometrics Toolbox. The syntax for the size function in the function "tprod" is also changed to make the code excutable on earlier versions of Matlab.

- [tRPCA](https://github.com/canyilu/Tensor-Robust-Principal-Component-Analysis-TRPCA) (Lu, Canyi, et al. "Tensor robust principal component analysis with a new tensor nuclear norm." IEEE transactions on pattern analysis and machine intelligence 42.4 (2019): 925-938.)
