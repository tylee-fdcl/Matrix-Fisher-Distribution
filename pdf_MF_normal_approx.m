function c = pdf_MF_normal_approx(s,type_approx,bool_scaled)
%pdf_MF_normal_approx: the approximated normalizing constant for the matrix Fisher distribution
%on SO(3)
%   c = pdf_MF_normal_approx(s) is the approximated normalizing constant for the 
%   matrix Fisher distribution on SO(3), for a given 3x1 (or 1x3) proper singular
%   values s.
%
%   c = pdf_MF_normal(s,TYPE_APPROX) returns the value 
%   specified by TYPE_APPROX:
%       0 - approximation by almost uniform distribuitons when s is small
%       1 - approximaiton by highly concentraed distributions when s_i+s_j
%       is large
%
%   c = pdf_MF_normal(s,TYPE_APPROX,BOOL_SCALED) returns the scaled value 
%   depending on BOOL_SCALED:
%       0 - (default) is the same as pdf_MF_normal(s,TYPE_APPROX)
%       1 - returnes an exponentially scaled normlaizing constant,
%       exp(-sum(s))*c
%
%   See T. Lee, "Bayesian Attitude Estimation with the Matrix Fisher
%   Distribution on SO(3)", 2017, http://arxiv.org/abs/1710.03746,
%   also T. Lee, "Bayesian Attitude Estimation with Approximate Matrix 
%   Fisher Distributions on SO(3)", 2018

assert(or(min(size(s)==[1 3]),min(size(s)==[3 1])),'ERROR: s should be 3 by 1 or 1 by 3');
assert(or(type_approx==1,type_approx==0),'ERROR: type_approx should be 0 or 1');

% if bool_scaled is not defined, then set it false
if nargin < 3
    bool_scaled=false;
end

switch type_approx
    case 0
        c=1+1/6*(s(1)^2+s(2)^2+s(3)^2)+1/6*s(1)*s(2)*s(3);
        if bool_scaled
            c=c/exp(sum(s));
        end
    case 1
        c=1/sqrt(8*pi*(s(1)+s(2))*(s(2)+s(3))*(s(3)+s(1)));
        if ~bool_scaled
            c=c*exp(sum(s));
        end
end
