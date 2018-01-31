function s=pdf_MF_M2S_approx(d,type_approx)
%pdf_MF_M2S_approx: transforms the first moments into the proper singular values
%   s=pdf_MF_M2S_approx(d,TYPE_APPROX) numerically solves the following equations for s
%
%       \frac{1}{c(S)}\frac{\partial c(S)}{s_i} - d_i = 0,
%
%   to find the proper singular values of the matrix Fisher distribution
%   whose first moments match with d. Using the approximated values for
%   c(S), the above equation is solved explicitly. 
%
%   The variable TYPE_APPROX determines the type of approximation:
%       0 - approximation by almost uniform distribuitons when s is small
%       1 - approximaiton by highly concentraed distributions when s_i+s_j
%       is large
%
%   See T. Lee, "Bayesian Attitude Estimation with the Matrix Fisher
%   Distribution on SO(3)", 2017, http://arxiv.org/abs/1710.03746,
%   also T. Lee, "Bayesian Attitude Estimation with Approximate Matrix
%   Fisher Distributions on SO(3)", 2018
%
%   See also PDF_MF_M2S

assert(or(min(size(d)==[1 3]),min(size(d)==[3 1])),'ERROR: d should be 3 by 1 or 1 by 3');
assert(or(type_approx==1,type_approx==0),'ERROR: type_approx should be 0 or 1');

switch type_approx
    case 0
        s=3*d;
    case 1
        s=zeros(3,1);
        for i=1:3
            index=circshift([1 2 3],[0 4-i]);
            j=index(2);
            k=index(3);
            
            s(i)=1/2*(-1/(1+d(i)-d(j)-d(k))+1/(1-d(i)+d(j)-d(k))+1/(1-d(i)-d(j)+d(k)));
        end
end
