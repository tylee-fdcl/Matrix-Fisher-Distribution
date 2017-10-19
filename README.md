# Matrix-Fisher-Distribution

The matrix Fisher distribution is a compact form of an exponential density model developed for random matrices [[1]](#Mar). This repository contains Matlab files to perform various stochastic analyses for the matrix Fisher distribution on the special orthogonal group SO(3). 

The mathematical fomulation of the presented algorithms are available at the following paper:

- T. Lee, ["*Bayesian Attitude Estimation with the Matrix Fisher Distribution on SO(3)*"](https://arxiv.org/abs/1710.03746/) 	arXiv:1710.03746, 2017


```pdf_MF.m : compute the probability density
pdf_MF_normal.m : compute the normalizing constant
pdf_MF_normal_deriv.m : compute the derivatives of the normalizing constant

pdf_MF_moment.m : compute the moment
pdf_MF_M2S.m : convert the first moment into the proper singular values

pdf_MF_unscented_transform.m : perform the unscented transform
pdf_MF_inv_unscented_transform.m : perform the inverse unscented transform
psvd.m : perform proper singular value decomposition
```

## Author
 - Taeyoung Lee, [Flight Dynamics and Control Lab](http://fdcl.seas.gwu.edu/), George Washington University 

## References
1.  <a name="Mar">K. Mardia and P. Jupp, *Directional Statistics,* Wiley, 1999.</a>

