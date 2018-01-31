function dc_out=pdf_MF_normal_deriv_approx(s,type_approx,bool_scaled)
%pdf_MF_norma_deriv_approx: the approximated derivatives of the normalizing constant for the matrix Fisher distribution
%on SO(3)
%   [dc] = pdf_MF_normal(s,TYPE_APPROX,BOOL_SCALED) returns the 3x1 first
%   order derivative dc of the normalizing constant with respect to the
%   proper singular values for the matrix Fisher distribution on SO(3),
%   for a given 3x1 (or 1x3) proper singular values s.
%
%   dc = pdf_MF_normal_deriv_approx(s,TYPE_APPROX) returns the value
%   specified by TYPE_APPROX:
%       0 - approximation by almost uniform distribuitons when s is small
%       1 - approximaiton by highly concentraed distributions when s_i+s_j
%       is large
%
%   c = pdf_MF_normal_deriv_approx(s,TYPE_APPROX,BOOL_SCALED) returns the scaled value
%   depending on BOOL_SCALED:
%       0 - (default) is the same as pdf_MF_normal_deriv_approx(s,TYPE_APPROX)
%       1 - computes the derivatives of the exponentially scaled normalizing constant,
%       c_bar = exp(-sum(s))*c
%
%   See T. Lee, "Bayesian Attitude Estimation with the Matrix Fisher
%   Distribution on SO(3)", 2017, http://arxiv.org/abs/1710.03746,
%   also T. Lee, "Bayesian Attitude Estimation with Approximate Matrix
%   Fisher Distributions on SO(3)", 2018
%
%   See also PDF_MF_NORMAL_DERIV

assert(or(min(size(s)==[1 3]),min(size(s)==[3 1])),'ERROR: s should be 3 by 1 or 1 by 3');
assert(or(type_approx==1,type_approx==0),'ERROR: type_approx should be 0 or 1');

% if bool_scaled is not defined, then set it false
if nargin < 2
    bool_scaled=false;
end


if ~bool_scaled
    c=pdf_MF_normal_approx(s,type_approx,bool_scaled);
    switch type_approx
        case 0
            dc=1/3*[s(1);s(2);s(3)]+1/6*[s(2)*s(3); s(3)*s(1); s(1)*s(2)];
        case 1
            dc = zeros(3,1);
            for i=1:3
                index=circshift([1 2 3],[0 4-i]);
                j=index(2);
                k=index(3);
                
                dc(i)=c*(1-0.5*(1/(s(i)+s(j))+1/(s(i)+s(k))));
            end
    end
    dc_out=dc;
else
    c_bar=pdf_MF_normal_approx(s,type_approx,bool_scaled);
    
    switch type_approx
        case 0            
            dc_bar=1/3*[s(1);s(2);s(3)]+1/6*[s(2)*s(3); s(3)*s(1); s(1)*s(2)]*exp(-sum(s))-c_bar;
        case 1
            dc_bar = zeros(3,1);
            
            for i=1:3
                index=circshift([1 2 3],[0 4-i]);
                j=index(2);
                k=index(3);
                
                dc_bar(i)=c_bar*(1-0.5*(1/(s(i)+s(j))+1/(s(i)+s(k))));
            end
    end
    
    dc_out=dc_bar;
    
end

