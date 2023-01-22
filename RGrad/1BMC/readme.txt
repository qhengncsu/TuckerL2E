One-Bit Matrix Completion (1BMC) Toolbox version 1.0

This software package contains tools for matrix completion in the extreme case of noisy 1-bit observations. For further details, see the paper "1-Bit Matrix Completion," by Mark A. Davenport, Yaniv Plan, Ewout van den Berg, and Mary Wootters at http://arxiv.org/abs/1209.3672.

This toolbox contains most of the software necessary to reproduce the results presented in this paper. In particular, see the file "syntheticSims.m" reproduces the synthetic simulations and "movielens.m" reproduces the simulations involving the MovieLens dataset. (Note that "movielens.m" also requires the installation of the TFOCS package and the raw MovieLens dataset, which are not included in this package.  See "movielens.m" for further details.)

To begin experimenting with one-bit matrix completion, see the "demo.m" file, which shows the core routines of this toolbox in action.  

Please e-mail mdav@gatech.edu if you find any bugs or have any questions. 