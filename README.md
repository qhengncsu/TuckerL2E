# Code for  "Robust Low-rank Tensor Decomposition with the L2 Criterion"

Qiang Heng, Eric C. Chi, Yufeng Liu (2023). [arXiv:2208.11806](https://arxiv.org/pdf/2208.11806.pdf)

## Dependencies
This repo is dependent on [Tensor Toolbox for MATLAB](https://www.tensortoolbox.org/) and the [L-BFGS-B implementation by Stephen Becker](https://github.com/stephenbeckr/L-BFGS-B-C).

Steps to get the code up and running:
1. Download and install [Matlab 2019b](https://www.mathworks.com/products/matlab.html) or newer, assuming that you currently do not have Matlab on your machine.
2. Install [Statistics and Machine Learning Toolbox](https://www.mathworks.com/products/statistics.html).
3. 

1. See demo_simulated_tensor.m for examples with simulated data, this corresponds to section 5.1 in the paper.
2. See demo_fMRI.m for the tensor denoising application, this corresponds to section 6.1 in the paper.
3. See demo_dorrit.m for the Chemometrics application, this corresponds to section 6.2 in the paper.

demo_dorrit.m will run relatively fast, demo_simulated_tensor.m and demo_fMRI.m will take about 5-10 minutes.
