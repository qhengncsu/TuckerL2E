# Code for  "Robust Low-rank Tensor Decomposition with the L2 Criterion"

Qiang Heng, Eric C. Chi, Yufeng Liu (2023).
This repository is maintained by Qiang Heng (qheng@ncsu.edu).

## Dependencies
- [Matlab 2019b](https://www.mathworks.com/products/matlab.html) or newer
- [Statistics and Machine Learning Toolbox](https://www.mathworks.com/products/statistics.html), which needs to be installed with the Add-on Explorer within Matlab.

Steps to get the code up and running:
1. Install [Tensor Toolbox for MATLAB v3.2.1](https://gitlab.com/tensors/tensor_toolbox/-/releases/v3.2.1). Please download the zip file [Source code (zip)](https://gitlab.com/tensors/tensor_toolbox/-/archive/v3.2.1/tensor_toolbox-v3.2.1.zip) and unzip it. Then please place the folder at "[path of repo]/TuckerL2E/tensor_toolbox-v3.2.1". An example absolute path on my machine is "C:/Users/Qiang/Documents/code_TuckerL2E/TuckerL2E/tensor_toolbox-v3.2.1", where "C:/Users/Qiang/Documents/code_TuckerL2E" is the path for this repo. 
2. Download as zip the github repo [L-BFGS-B-C](https://github.com/stephenbeckr/L-BFGS-B-C). Unzip the zip file. Notice that the folder name will be "L-BFGS-B-C-master" after being downloaded from github. Please rename the folder to "L-BFGS-B-C" and place the folder at "[path of repo]/TuckerL2E/L-BFGS-B-C". An example absolute path on my machine is "C:/Users/Qiang/Documents/code_TuckerL2E/TuckerL2E/L-BFGS-B-C". When you are done, the organization of files should look like [this](https://github.com/qhengncsu/TuckerL2E/blob/master/images/layout.png) within Matlab.
3. Please rename the "hosvd.m" function inside "[path of repo]/TuckerL2E/tensor_toolbox-v3.2.1" to "hosvd_tt.m". See [here](https://github.com/qhengncsu/TuckerL2E/blob/master/images/hosvd_tt.png) for a preview. The function name written inside "hosvd_tt.m" must be changed as well. In other words, the first line of "hosvd_tt.m" should now be "function T = hosvd_tt(X,tol,varargin)". 


## Compiling the mex files
The [L-BFGS-B-C](https://github.com/stephenbeckr/L-BFGS-B-C) repo requires compilation! Thus you must follow the following steps to compile the mex files:

4. For Windows users, please install "MATLAB Support for MinGW-w64 C/C++ Compiler" using the Add-on Explorer within Matlab.
5. For Mac users, please ensure Xcode is installed and functional.
6. For Linux users, make sure GCC compiler is installed (most linux distributions have GCC pre-installed).
7. Navigate to "[path of repo]/TuckerL2E/L-BFGS-B-C/Matlab" and run the script "compile_mex.m" within Matlab (The working directory of Matlab should be "[path of repo]/TuckerL2E/L-BFGS-B-C/Matlab"). 

## Example scripts
Now you are ready to try out the scripts! (Remember to swith the working directory back to [path of repo] after step 7)

8. "generate_figure_2_subplot_[i].m" (i=1,2,3,4) will help you recreate the four subplots of Figure 2 in the paper. 
9. "generate_figure_6.m" will help you recreate Figure 6 in the paper. 
10. "generate_figure_8.m" will help you recreate Figure 8 in the paper. 

"generate_figure_8.m" will run relatively fast, "generate_figure_6.m" will take about 10-15 minutes, each one of "generate_figure_2_subplot_[i].m" will take about 15-25 minutes if the number of replicate "N" is set to 1. Notice that in the paper, Figure 2 is generated with 50 random replicates. However, one random replicate is enough to generate a figure that appears similar to Figure 2. 
