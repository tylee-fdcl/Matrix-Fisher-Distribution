function [F s]=pdf_MF_inv_unscented_transform(R,W0,W)
%pdf_MF_inv_ unscented_transform: inverse unscented transform for the matrix 
% Fisher distribution on SO(3)
%   [F s]=pdf_MF_inv_unscented_transform(R,W0,W) performs the inverse 
%   unscented transform for the matrix Fisher distribution on SO(3) to
%   convert the sigma points and weights to a matrix parameter whose first
%   moment is equal to the weighted sum of the sigma points. 
%   It also returns the proper singular value of F.
%
%   The input variables are defined as
%       R (3x3x7) - sigma points ordered as R0, R1+, R1-, R2+, R2-, R3+, R3-
%       W0 (1x1) - weighting for R0
%       W (3x1) - weighting for R1+-, R2+-, R3+-
%
%   See also PDF_MF_UNSCENTED_TRANSFORM.
%
%   See T. Lee, "Bayesian Attitude Estimation with the Matrix Fisher
%   Distribution on SO(3)", 2017, http://arxiv.org/abs/1710.03746

M1=W0*R(:,:,1);
for i=1:3
    M1=M1+W(i)*R(:,:,2*i)+W(i)*R(:,:,2*i+1);
end

[U D V]=psvd(M1);

s=pdf_MF_M2S(diag(D));

F=U*diag(s)*V';
