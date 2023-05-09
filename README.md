# Code for  "Robust Low-rank Tensor Decomposition with the L2 Criterion"

Qiang Heng, Eric C. Chi, Yufeng Liu (2023). doi: [10.1080/00401706.2023.2200541](https://www.tandfonline.com/doi/abs/10.1080/00401706.2023.2200541)

This repository is maintained by Qiang Heng (qheng@ncsu.edu).


## Dependencies
- [Matlab 2019b](https://www.mathworks.com/products/matlab.html) or newer
- [Statistics and Machine Learning Toolbox](https://www.mathworks.com/products/statistics.html), which needs to be installed with the Add-on Explorer within Matlab.

## Compiling the mex files
The [L-BFGS-B-C](https://github.com/stephenbeckr/L-BFGS-B-C) repo requires compilation! Thus you must follow the following steps to compile the mex files:

1. For Windows users, please install "MATLAB Support for MinGW-w64 C/C++ Compiler" using the Add-on Explorer within Matlab.
2. For Mac users, please ensure Xcode is installed and functional.
3. For Linux users, make sure GCC compiler is installed (most linux distributions have GCC pre-installed).
4. **Navigate to "[path of repo]/TuckerL2E/L-BFGS-B-C/Matlab" and run the script "compile_mex.m" within Matlab (The working directory of Matlab should be "[path of repo]/TuckerL2E/L-BFGS-B-C/Matlab").** 

## Example scripts
Now you are ready to try out the scripts! (Remember to switch the working directory back to [path of repo] after step 4)

5. "generate_figure_2_subplot_[i].m" (i=1,2,3,4) will help you recreate the four subplots of Figure 2 in the paper. 
6. "generate_figure_6.m" will help you recreate Figure 6 in the paper. 
7. "generate_figure_8.m" will help you recreate Figure 8 in the paper. 

"generate_figure_8.m" will run relatively fast, "generate_figure_6.m" will take about 10-15 minutes, each one of "generate_figure_2_subplot_[i].m" will take about 15-25 minutes if the number of replicate "N" is set to 1. Notice that in the paper, Figure 2 is generated with 50 random replicates. However, one random replicate is enough to generate a figure that appears similar to Figure 2. 
