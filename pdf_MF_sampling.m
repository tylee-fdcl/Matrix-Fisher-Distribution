function [R accept_ratio]=pdf_MF_sampling(F,N)
%pdf_MF_sampling: samping for the matrix Fisher distribution on SO(3)
%   R=pdf_MF_sampling(F,N) returns N rotation matricies distributed
%   according to the matrix Fisher distribution with hte matrix parameter F
%
%   See T. Lee, "Bayesian Attitude Estimation with the Matrix Fisher
%   Distribution on SO(3)", 2017, http://arxiv.org/abs/1710.03746

[U S V]=psvd(F);
B=zeros(4,4);
B(1:3,1:3)=2*S-trace(S)*eye(3);
B(4,4)=trace(S);

[x accept_ratio]=sampBing(B,N);

R=zeros(3,3,N);
for i=1:N
    q4=x(4,i);
    q=x(1:3,i);
    R(:,:,i)=(q4^2-q'*q)*eye(3)+2*q*q'+2*q4*hat(q);
    R(:,:,i)=U*R(:,:,i)*V';
end

end

function [x accept_ratio]=sampBing(B,N)
% simulating Bingham distribution
% See Kent, Ganeiber, and Mardia, "A new method to simulate the Bingham and
% related distribution, 2013

lamB=eig(B);
A=-B;
lamA=-lamB;
min_lamA=min(lamA);
lamA=lamA-min_lamA;
A=A-min_lamA*eye(4);

funb = @(b) 1/(b+2*lamA(1))+1/(b+2*lamA(2))+1/(b+2*lamA(3))+1/(b+2*lamA(4))-1;

tol = optimoptions('fsolve','TolFun', 1e-8, 'TolX', 1e-8,'display','off');
[b err exitflag]=fsolve(funb,1,tol);
if exitflag ~= 1
    disp([err exitflag]);
end

W=eye(4)+2*A/b;
Mstar=exp(-(4-b)/2)*(4/b)^2;

x=zeros(4,N);
nx=0;
nxi=0;
while nx < N
    xi=mvnrnd(zeros(4,1),inv(W))';
    xi=xi/norm(xi);
    nxi=nxi+1;
    
    pstar_Bing=exp(-xi'*A*xi);
    pstar_ACGD=(xi'*W*xi)^(-2);
    u=rand(1);
    
    if u < (pstar_Bing / (Mstar*pstar_ACGD))        
        nx=nx+1;
        x(:,nx)=xi;        
    end
    
end

accept_ratio=N/nxi;
end
