LASSO-Solver-OMP
=========================

>**Designer:** Junbo Zhao, Wuhan University, Working in Tsinghua National lab of intelligent images and documents processing.      
**Contact:** zhaojunbo1992chasing@gmail.com +86-18672365683     

**Introduction**     
This package implements Orthonogal Matching Pursuit algorithm (OMP) approach as a famous LASSO solver, and this program is written in C++ assisted by LAPACK.     
LASSO is a critical issue, which can be treated as a statistic problem but is widely exploited in lots of applications. Sparse coding, for example as a important tool for computer vision, natural language processing and machine learning, is based on the well-development of LASSO solvers.
OMP is well-known, due to its advantages over Basis pursuit or previously proposed method Matching Pursuit (MP). OMP achieves a quicker convergence, and overcomes some drawbacks of the other methods.
Many codes, implementing OMP algorithm, are mostly implemented in MATLAB, for its convenience. However, since MATLAB is not so advantageous when encountering big iterations, for some large scale or high dimensional problems, C++ is prefered.
The specific introduction of this OMP algorithm could be seen in the paper *Orthogonal Matching Pursuit-Recursive Function Approximation with Applications to wavelet decomposition*, 2003.     

**Configuration**    
To run the project, you should pre-configure LAPACK interfaces in your Visual C++ environment. It will be quite easy if you follow the website: http://icl.eecs.utk.edu/lapack-for-windows/clapack/index.html#libraries      
Find the "Part2: Using CLAPACK subroutines in a Visual (Studio) C/C++ Project" section in chapter "Running CLAPACK under Windows", and please follow the instructions strictly. Moreover, files will not be included in the repository if they can be derived from the former website (namely some .dll and .lib files).    
**Also, I will not upload some header files for CLAPACK in this repository and you can download them in LDA repository, which I constructed some days earlier.**     

**Platform**     
This program is tested on VS2010, 32bit, Win7 system. I cannot guarantee it could be adopted on other platforms. If you can successfully compile the source code on some other platforms, please feel free to contact me！     

**Usage**     
1. The mechanism of the matrix operations and storage is quite similar to previous repositories. You can have a look at the readme.md file of PCA or LDA repositories if you need some explanations.    
2. Specific usages and the meaning of parameters can be found in omp.h
