function varargout=pdf_MF_normal_deriv(s,bool_ddc,bool_scaled)
%pdf_MF_norma_deriv: the derivatives of the normalizing constant for the matrix Fisher distribution
%on SO(3)
%   [dc, ddc] = pdf_MF_normal(s,BOOL_DDC,BOOL_SCALED) returns the 3x1 first
%   order derivative dc and the 3x3 second order derivatives ddc of the normalizing
%   constant with respect to the proper singular values for the matrix Fisher
%   distribution on SO(3), for a given 3x1 (or 1x3) proper singular
%   values s.
%
%   BOOL_DDC determines whether the second order derivative
%   are computed or not:
%       0 - (defalut) is the same as dc=pdf_MF_normal_deriv(s), and the second
%       order derivatives are not computed
%       1 - computes the second order derivatives, and reuturns ddc
%
%   BOOL_SCALED determines whether the normalizing constant is
%   scaled or not:
%       0 - (defalut) is the same as pdf_MF_normal_deriv(s,BOOL_DDC), and
%       then derivatives of the unscaled normalizing constant c are returned
%       1 - computes the derivatives of the exponentially scaled normalizing constant,
%       c_bar = exp(-sum(s))*c
%
%   Examples
%       dc=pdf_MF_normal_deriv(s) - first order derivatives of c
%       [dc, ddc]=pdf_MF_normal_deriv(s,true) - first and second
%       order derivatives of c
%
%       dc_bar=pdf_MF_normal_deriv(s,false,true) - first order
%       derivatives of the exponentially scaled c
%       [dc_bar, ddc_bar]=pdf_MF_normal_deriv(s,true,true) - first and second
%       order derivatives of the exponentially scaled c
%
%   See T. Lee, "Bayesian Attitude Estimation with the Matrix Fisher
%   Distribution on SO(3)", 2017, http://arxiv.org/abs/1710.03746


if nargin < 3
    bool_scaled = false;
end
if nargin < 2
    bool_ddc = false;
end

if ~bool_scaled
    dc=zeros(3,1);
    
    % derivatives of the normalizing constant
    for i=1:3
        dc(i) = integral(@(u) f_kunze_s_deriv_i(u,s,i),-1,1);
    end
    varargout{1}=dc;
    
    if bool_ddc
        % compute the second order derivatives of the normalizing constant
        ddc=zeros(3,3);
        for i=1:3
            ddc(i,i) = integral(@(u) f_kunze_s_deriv_ii(u,s,i),-1,1);
            for j=i+1:3
                ddc(i,j) = integral(@(u) f_kunze_s_deriv_ij(u,s,i,j),-1,1);
                ddc(j,i)=ddc(i,j);
            end
        end
        
        varargout{2}=ddc;
    end
    
else
    % derivatives of the scaled normalizing constant
    dc_bar = zeros(3,1);
    
    for i=1:3
        index=circshift([1 2 3],[0 4-i]);
        j=index(2);
        k=index(3);
        
        dc_bar(k) = integral(@(u) f_kunze_s_deriv_scaled(u,[s(i),s(j),s(k)]),-1,1);
    end
    varargout{1}=dc_bar;
    
    % compute the second order derivatives of the scaled normalizing
    % constant
    if bool_ddc
        ddc_bar=zeros(3,3);
        for i=1:3
            index=circshift([1 2 3],[0 4-i]);
            j=index(2);
            k=index(3);
            
            c_bar=pdf_MF_normal(s,1);
            ddc_bar(k,k) = integral(@(u) f_kunze_s_deriv_scaled_kk(u,[s(i),s(j),s(k)]),-1,1);
            ddc_bar(i,j) = integral(@(u) f_kunze_s_deriv_scaled_ij(u,[s(i),s(j),s(k)]),-1,1) ...
                -dc_bar(i)-dc_bar(j)-c_bar;
            ddc_bar(j,i) = ddc_bar(i,j);
        end
        varargout{2}=ddc_bar;
    end
end
end

function Y=f_kunze_s_deriv_scaled(u,s)
% integrand for the derivative of the scaled normalizing constant
[l m]=size(u);
for ii=1:l
    for jj=1:m
        J=besseli(0,1/2*(s(1)-s(2))*(1-u(ii,jj)),1)*besseli(0,1/2*(s(1)+s(2))*(1+u(ii,jj)),1);
        Y(ii,jj)=1/2*J*(u(ii,jj)-1)*exp((min([s(1) s(2)])+s(3))*(u(ii,jj)-1));
    end
end
end


function Y=f_kunze_s_deriv_scaled_kk(u,s)
% integrand for the second order derivative of the scaled normalizing constant
[l m]=size(u);
for ii=1:l
    for jj=1:m
        J=besseli(0,1/2*(s(1)-s(2))*(1-u(ii,jj)),1)*besseli(0,1/2*(s(1)+s(2))*(1+u(ii,jj)),1);
        Y(ii,jj)=1/2*J*(u(ii,jj)-1)^2*exp((min([s(1) s(2)])+s(3))*(u(ii,jj)-1));
    end
end
end

function Y=f_kunze_s_deriv_scaled_ij(u,s)
% integrand for the scaled second order derivative of the normalizing constant
[l m]=size(u);
for ii=1:l
    for jj=1:m
        J=1/4*besseli(1,1/2*(s(2)-s(3))*(1-u(ii,jj)),1)*besseli(0,1/2*(s(2)+s(3))*(1+u(ii,jj)),1) ...
            *exp((s(1)+min([s(2) s(3)]))*(u(ii,jj)-1))*u(ii,jj)*(1-u(ii,jj)) ...
            +1/4*besseli(0,1/2*(s(2)-s(3))*(1-u(ii,jj)),1)*besseli(1,1/2*(s(2)+s(3))*(1+u(ii,jj)),1) ...
            *exp((s(1)+min([s(2) s(3)]))*(u(ii,jj)-1))*u(ii,jj)*(1+u(ii,jj));
        Y(ii,jj)=J;
    end
end
end

function Y=f_kunze_s_deriv_i(u,s,i)
% integrand for the derivative of the normalizing constant
index=circshift([1 2 3],[0 4-i]);
j=index(2);
k=index(3);

[l m]=size(u);
for ii=1:l
    for jj=1:m
        J00=besseli(0,1/2*(s(j)-s(k))*(1-u(ii,jj)))*besseli(0,1/2*(s(j)+s(k))*(1+u(ii,jj)));
        Y(ii,jj)=1/2*J00*u(ii,jj)*exp(s(i)*u(ii,jj));
    end
end
end

function Y=f_kunze_s_deriv_ii(u,s,i)
% integrand for the second-order derivative of the normalizing constant
index=circshift([1 2 3],[0 4-i]);
j=index(2);
k=index(3);

[l m]=size(u);
for ii=1:l
    for jj=1:m
        J00=besseli(0,1/2*(s(j)-s(k))*(1-u(ii,jj)))*besseli(0,1/2*(s(j)+s(k))*(1+u(ii,jj)));
        Y(ii,jj)=1/2*J00*u(ii,jj)^2*exp(s(i)*u(ii,jj));
    end
end
end

function Y=f_kunze_s_deriv_ij(u,s,i,j)
% integrand for the mixed second-order derivative of the normalizing constant
k=setdiff([1 2 3],[i,j]);

[l m]=size(u);
for ii=1:l
    for jj=1:m
        J10=besseli(1,1/2*(s(j)-s(k))*(1-u(ii,jj)))*besseli(0,1/2*(s(j)+s(k))*(1+u(ii,jj)));
        J01=besseli(0,1/2*(s(j)-s(k))*(1-u(ii,jj)))*besseli(1,1/2*(s(j)+s(k))*(1+u(ii,jj)));
        Y(ii,jj)=1/4*J10*u(ii,jj)*(1-u(ii,jj))*exp(s(i)*u(ii,jj))...
            +1/4*J01*u(ii,jj)*(1+u(ii,jj))*exp(s(i)*u(ii,jj));
    end
end
end
