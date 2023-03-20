# Code for  "Robust Low-rank Tensor Decomposition with the L2 Criterion"

Qiang Heng, Eric C. Chi, Yufeng Liu (2023).
This repository is maintained by Qiang Heng (qheng@ncsu.edu).

## Dependencies
- [Matlab 2019b](https://www.mathworks.com/products/matlab.html) or newer
- [Statistics and Machine Learning Toolbox](https://www.mathworks.com/products/statistics.html)
This repo is also dependent on [Tensor Toolbox for MATLAB](https://www.tensortoolbox.org/) and the [L-BFGS-B implementation by Stephen Becker](https://github.com/stephenbeckr/L-BFGS-B-C).

## Comipling the mex files
After downloading and unzipping this repo, you will also need to compile the mex files for [L-BFGS-B implementation by Stephen Becker](https://github.com/stephenbeckr/L-BFGS-B-C) before you attempt to run the scripts.

1. For Windows users, please install "MATLAB Support for MinGW-w64 C/C++ Compiler" using the Add-on Explorer within Matlab.
2. For Mac users, please ensure Xcode is installed and functional.
3. For Linux users, the binary files included should in general be directly usable. If not, go to step 4.
4. Navigate to "./TuckerL2E/L-BFGS-B-C/Matlab" and run the script "compile_mex.m" within Matlab (The working directory of Matlab should be "./TuckerL2E/L-BFGS-B-C/Matlab"). 

## Example Scripts
Now you are ready to try out the scripts!

1. "generate_figure_2.m" this will help you recreate the fourth subplot of Figure 2 in the paper (using only one random replicate). 
2. "generate_figure_6.m" this will help you recreate Figure 6 in the paper. 
3. "generate_figure_8.m" this will help you recreate Figure 8 in the paper. 

"generate_figure_8.m" will run relatively fast, "generate_figure_6.m" will take about 5-10 minutes, "generate_figure_2.m" will take about 15-25 minutes.
