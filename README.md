# Code for  "Robust Low-rank Tensor Decomposition with the L2 Criterion"

Qiang Heng, Eric C. Chi, Yufeng Liu (2023). [arXiv:2208.11806](https://arxiv.org/pdf/2208.11806.pdf)

## Dependencies
This repo is dependent on [Tensor Toolbox for MATLAB](https://www.tensortoolbox.org/) and the [L-BFGS-B implementation by Stephen Becker](https://github.com/stephenbeckr/L-BFGS-B-C).

Steps to get the code up and running (step 1 and 2 can be skipped if you already have them):
1. Download and install [Matlab 2019b](https://www.mathworks.com/products/matlab.html) or newer.
2. Install [Statistics and Machine Learning Toolbox](https://www.mathworks.com/products/statistics.html).
3. Download this repository and unzip it to somewhere you like. (If you downloaded the repo as zip on github, the folder name will be "TuckerL2E-master", we first have to rename the folder to "TuckerL2E")
4. Download [tensor_toolbox-v3.2.1](https://gitlab.com/tensors/tensor_toolbox/-/releases/v3.2.1), unzip the file and place the folder in this repository such that we can access its contents with the relative path "./TuckerL2E/tensor_toolbox-v3.2.1". An example absolute path would be "C:/Users/Qiang/Desktop/TuckerL2E/TuckerL2E/tensor_toolbox-v3.2.1".
5. Download [L-BFGS-B implementation by Stephen Becker](https://github.com/stephenbeckr/L-BFGS-B-C), unzip the file and place the folder in this repository such that we can access its contents with the relative path "./TuckerL2E/L-BFGS-B-C" (Again, if you download as zip on github, the folder name will be "L-BFGS-B-C-master", we first have to rename the folder to lose the "-master").
6. Navigate to "./TuckerL2E/L-BFGS-B-C/Matlab" and run the script "compile_mex.m" within Matlab. This will you help compile the mex files provided.

Now you are ready to try out the scripts!

1. See demo_simulated_tensor.m for examples with simulated data, this corresponds to section 5.1 in the paper.
2. See demo_fMRI.m for the tensor denoising application, this corresponds to section 6.1 in the paper.
3. See demo_dorrit.m for the Chemometrics application, this corresponds to section 6.2 in the paper.

demo_dorrit.m will run relatively fast, demo_simulated_tensor.m and demo_fMRI.m will take about 5-10 minutes.
