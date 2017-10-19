function p=pdf_MF(R,F)
%pdf_MF: the probability density for the matrix Fisher distribution
%on SO(3)
%   p = pdf_MF(R,F) is the probabilty density for the matrix Fisher
%   distribution on SO(3) with the matrix parameter F at the attitude
%   represented by the rotation matrix R in SO(3).
%
%   See T. Lee, "Bayesian Attitude Estimation with the Matrix Fisher
%   Distribution on SO(3)", 2017, http://arxiv.org/abs/1710.03746

[U S V]=psvd(F);

% nominal implementation
% c=pdf_MF_normal(diag(S));
% p=1/c*exp(trace(F'*R))

% scaled implementation
c_bar=pdf_MF_normal(diag(S),1);
p=1/c_bar*exp(trace(F'*R)-trace(S));

end