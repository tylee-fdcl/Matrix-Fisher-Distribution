# Matrix-Fisher-Distribution

The matrix Fisher distribution is a compact form of an exponential density model developed for random matrices [[1]](#Mar). This repository contains Matlab files to perform various stochastic analyses and attitude estimation for the matrix Fisher distribution on the special orthogonal group SO(3). 

The mathematical fomulation of the presented algorithms are available at the following paper:

- T. Lee, ["*Bayesian Attitude Estimation with the Matrix Fisher Distribution on SO(3)*"](https://arxiv.org/abs/1710.03746/) 	arXiv:1710.03746, 2017

## List of Functions

```pdf_MF.m : compute the probability density
pdf_MF_normal.m : compute the normalizing constant
pdf_MF_normal_deriv.m : compute the derivatives of the normalizing constant

pdf_MF_moment.m : compute the moment
pdf_MF_M2S.m : convert the first moment into the proper singular values

pdf_MF_sampling.m : generate sample attitudes

pdf_MF_unscented_transform.m : perform the unscented transform
pdf_MF_inv_unscented_transform.m : perform the inverse unscented transform

psvd.m : perform proper singular value decomposition

test_pdf_MF.m : test the package

est_MF.m : attitude estimation
```
## Installing and Testing
Download all of the `m` files into a folder. Run `test_pdf_MF` in the Matlab workspace, and the results should be as follows.

```>> test_pdf_MF
check:c_bar    1.1543e-16
check:dc_bar   7.0903e-16
check:ddc_bar  8.9409e-16
check:M1       3.1032e-16
check:M2       8.0138e-16
check:M2S      1.2266e-15
check:UT       1.2542e-15
check:invUT    8.9461e-14
```

## Author
 - Taeyoung Lee, [Flight Dynamics and Control Lab](http://fdcl.seas.gwu.edu/), George Washington University 

## References
1.  <a name="Mar">K. Mardia and P. Jupp, *Directional Statistics,* Wiley, 1999.</a>

