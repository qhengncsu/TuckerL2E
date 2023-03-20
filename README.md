# Code for  "Robust Low-rank Tensor Decomposition with the L2 Criterion"

Qiang Heng, Eric C. Chi, Yufeng Liu (2023).
This repository is maintained by Qiang Heng (qheng@ncsu.edu).

## Dependencies
This repo is dependent on [Tensor Toolbox for MATLAB](https://www.tensortoolbox.org/) and the [L-BFGS-B implementation by Stephen Becker](https://github.com/stephenbeckr/L-BFGS-B-C).

Steps to get the code up and running (step 1 and 2 can be skipped if you already have them):
1. Download and install [Matlab 2019b](https://www.mathworks.com/products/matlab.html) or newer.
2. Install [Statistics and Machine Learning Toolbox](https://www.mathworks.com/products/statistics.html).
3. Download this repository and unzip it to somewhere you like (If you downloaded this repo as zip file from github, the folder name will be "TuckerL2E-master", we first have to rename the folder to "TuckerL2E").
4. Download [tensor_toolbox-v3.2.1](https://gitlab.com/tensors/tensor_toolbox/-/releases/v3.2.1), unzip the file and place the folder in this repository such that we can access its contents with the relative path "./TuckerL2E/tensor_toolbox-v3.2.1". An example absolute path would be "C:/Users/Qiang/Desktop/TuckerL2E/TuckerL2E/tensor_toolbox-v3.2.1".
5. Download [L-BFGS-B implementation by Stephen Becker](https://github.com/stephenbeckr/L-BFGS-B-C), unzip the file and place the folder in this repository such that we can access its contents with the relative path "./TuckerL2E/L-BFGS-B-C" (Again, if you download as zip from github, the folder name will be "L-BFGS-B-C-master", we first have to rename the folder to lose the "-master").

## Comipling the mex files
i. For Windows users, please install "MATLAB Support for MinGW-w64 C/C++ Compiler" using the Add-on Explorer within Matlab.
ii. For Mac users, please ensure you have Xcode installed.
iii. For Linux users, the binary files included should be directly usable.
iv. For Windows and Mac users, navigate to "./TuckerL2E/L-BFGS-B-C/Matlab" and run the script "compile_mex.m" within Matlab (The working directory of Matlab should be "./TuckerL2E/L-BFGS-B-C/Matlab"). 

## Example Scripts
Now you are ready to try out the scripts!

1. "generate_figure_2.m" this will help you recreate the fourth subplot of Figure 2 in the paper (using only one random replicate). 
2. "generate_figure_6.m" this will help you recreate Figure 6 in the paper. 
3. "generate_figure_8.m" this will help you recreate Figure 8 in the paper. 

"generate_figure_8.m" will run relatively fast, "generate_figure_6.m" will take about 5-10 minutes.
