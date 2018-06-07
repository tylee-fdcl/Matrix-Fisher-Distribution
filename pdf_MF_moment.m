function varargout=pdf_MF_moment(s,bool_M2)
%pdf_MF_moment: the canonical moment of the matrix Fisher distribution on SO(3)
%   [M1, M2, c_bar] = pdf_MF_normal(s,BOOL_M2) returns the 3x1 matrix M1 of 
%   the first moments and the 3x3 matrix of the non-zero second order moments 
%   for the matrix Fisher distribution with the parameter S=diag(s). 
%   It also returns the exponentially scaled normalizing constant. 
%
%   To obtain the 3x3 first moment for an arbitrary matrix parameter F, run 
%       [U S V]=psvd(F);
%       M1=U*diag(pdf_MF_moment(diag(S)))*V';
%
%   BOOL_M2 determines whether the second order canonical moments
%   are computed or not:
%       0 - (defalut) is the same as [M1, c_bar]=pdf_MF_moment(s), and the 
%       second order moments are not computed
%       1 - computes the second order moments
%
%   See T. Lee, "Bayesian Attitude Estimation with the Matrix Fisher
%   Distribution on SO(3)", 2017, http://arxiv.org/abs/1710.03746


if nargin < 2
    bool_M2 = false;
end

c_bar=pdf_MF_normal(s,1);

if ~bool_M2
    % compute the first order moments only
    
    dc_bar=pdf_MF_normal_deriv(s,0,1);
    M1=dc_bar/c_bar+1;

    varargout{1}=M1;
    varargout{2}=c_bar;
    
else
    % compute the first order moments and the second order moments
    
    [dc_bar, ddc_bar]=pdf_MF_normal_deriv(s,1,1);
    M1=dc_bar/c_bar+1;
    
    for i=1:3
        for j=i:3
            M2(i,j)=1+1/c_bar*(dc_bar(i)+dc_bar(j)+ddc_bar(i,j));
            M2(j,i)=M2(i,j);
        end
    end
    varargout{1}=M1;
    varargout{2}=M2;
    varargout{3}=c_bar;
end

end
