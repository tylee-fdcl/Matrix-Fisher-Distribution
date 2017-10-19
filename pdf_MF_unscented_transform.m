function [R W0 W theta sigma_min sigma]=pdf_MF_unscented_transform(F,sigma)
%pdf_MF_unscented_transform: unscented transform for the matrix Fisher
%distribution on SO(3)
%   [R W0 W theta sigma_min sigma]=pdf_MF_unscented_transform(F,sigma) performs the
%   unscented transform for the matrix Fisher distribution on SO(3) with
%   the 3x3 matrix parameter F with a scalar sigma satisfying the
%   following inequality,
%
%       sigma_min < sigma < 1,
%
%   where sigma_min = =max([(2*s(1)+s(2)-s(3)-1)/(2*s(1)+s(2)-s(3)+1),...
%   (s(1)-s(3))/(s(1)+s(2))]), and s is the proper sinagular value of F.
%
%   [R W0 W]=pdf_MF_unscented_transform(F) is equal to
%   [R W0 W]=pdf_MF_unscented_transform(F,(1+sigma_min)/2), which satisfies 
%   the above inqualiy for sigma by definition.
%
%   The function returns
%
%       R (3x3x7) - sigma points ordered as R0, R1+, R1-, R2+, R2-, R3+, R3-
%       W0 (1x1) - weighting for R0
%       W (3x1) - weighting for R1+-, R2+-, R3+-
%       theta (3x1) - (optional) rotation angles for R1+-, R2+-, R3+-
%       sigma_min (1x1) - (optional) lower bound of sigma
%       sigma (1x1) - (optional) sigma 
%
%   so that the weighted sum of the sigma points is equal to the first
%   moment of the matrix Fisher distribution.
%
%   See also PDF_MF_INV_UNSCENTED_TRANSFORM.
%
%   See T. Lee, "Bayesian Attitude Estimation with the Matrix Fisher
%   Distribution on SO(3)", 2017, http://arxiv.org/abs/1710.03746


[U S V]=psvd(F);
M=U*V';
s=diag(S);;
[M1 c_bar]=pdf_MF_moment(s);

sigma_min=max([(2*s(1)+s(2)-s(3)-1)/(2*s(1)+s(2)-s(3)+1),...
    (s(1)-s(3))/(s(1)+s(2))]);

if nargin < 2
    sigma=(sigma_min+1)/2;
else
    assert(and(sigma<1,sigma_min<sigma),...
        ['ERROR:the inequality for sigma is not satisfied. ' num2str(sigma_min,4) '<sigma<1, but sigma=' num2str(sigma,4)]);
end

E=eye(3);
R=zeros(3,3,7);
W=zeros(3,1);
R(:,:,1)=M;
for i=1:3
    index=circshift([1 2 3],[0 4-i]);
    j=index(2);
    k=index(3);
    
    if s(j)+s(k) >= 1
        theta(i)=acos(sigma+(1-sigma)*(log(c_bar)+sum(s)-s(i))/(s(j)+s(k)));
    else
        theta(i)=acos((sigma+(1-sigma)*(log(c_bar)+sum(s)-s(i))+1/2)*(s(j)+s(k))-1/2);
    end
    W(i)=1/4/(1-cos(theta(i)))*(M1(i)-M1(j)-M1(k)+1);
    R(:,:,2*i)=U*expm(theta(i)*hat(E(:,i)))*V';
    R(:,:,2*i+1)=U*expm(-theta(i)*hat(E(:,i)))*V';
end
W0=1-2*sum(W);

end
