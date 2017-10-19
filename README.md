# Matrix-Fisher-Distribution

The matrix Fisher distribution is a compact form of an exponential density model developed for random matrices [[1]](#Mar). This repository contains Matlab files to perform various stochastic analyses for the matrix Fisher distribution on the special orthogonal group SO(3). 

The mathematical fomulation of the presented algorithms are available at the following paper:

- T. Lee, ["*Bayesian Attitude Estimation with the Matrix Fisher Distribution on SO(3)*"](https://arxiv.org/abs/1710.03746/) 	arXiv:1710.03746, 2017


`pdf_MF.m
pdf_MF_M2S.m
pdf_MF_deriv.m
pdf_MF_inv_unscented_transform.m
pdf_MF_moment.m
pdf_MF_normal.m
pdf_MF_normal_deriv.m
pdf_MF_unscented_transform.m
psvd.m`

## References
1.  <a name="Mar">K. Mardia and P. Jupp, *Directional Statistics,* Wiley, 1999.</a>

