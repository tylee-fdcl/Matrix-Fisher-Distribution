function est_MF
%est_MF: attitude estimation with the matrix Fisher
%distributionp on SO(3)
%
%   Internal variables
%       - filename   : the name of the mat file where estimation results are saved
%       - EST_METHOD : determine the estimation scheme
%               0    : first order estimation
%               1    : unscented estimation
%       - F          : estimated matrix paramter
%       - s          : estimated proper singular values
%       - RM         : estimated mean attitude 
%       - rot_est_err : estimtion error in degrees
%
%   See T. Lee, "Bayesian Attitude Estimation with the Matrix Fisher
%   Distribution on SO(3)", 2017, http://arxiv.org/abs/1710.03746
close all;
global h J Jd m g e1 e2 e3 rho 

%% Estimation method and the initial estimate
filename='est_MF_first_order_0';
EST_METHOD=0;
F0=100*expm(pi*hat(e1));

% filename='est_MF_first_order_1';
% EST_METHOD=0;
% F0=zeros(3,3);
% % 
% filename='est_MF_unscented_0';
% EST_METHOD=1;
% F0=100*expm(pi*hat(e1));
% % 
% filename='est_MF_unscented_1';
% EST_METHOD=1;
% F0=zeros(3,3);

assert(EST_METHOD==0 || EST_METHOD==1,'EST_METHOD should be either 0 or 1');


if EST_METHOD==0
    disp('ESTIMATOR: First order attitude estimation');
else
    disp('ESTIMATOR: Unscented attitude estimation');
end

%% Rigid body parameters
m=1;
g=9.81;
e1=[1 0 0]';
e2=[0 1 0]';
e3=[0 0 1]';

Body.a=0.8;
Body.b=0.2;
Body.h=0.6;

rho=0.5*Body.h*e3;
J=diag([1/4*Body.b^2+1/3*Body.h^2,...
    1/4*Body.a^2+1/3*Body.h^2,...
    1/4*(Body.a^2+Body.b^2)]);
Jd=trace(J)/2*eye(3)-J;

%% True Trajectory
disp('Generating the true reference trajectories...');
h=0.02;
T=10;
t=0:h:T;
N=length(t);

R_true=zeros(3,3,N);
W_true=zeros(3,N);

R_true(:,:,1)=eye(3);
W_true(:,1)=1.3*sqrt(m*g*norm(rho)*2/trace(J))*[1 1 1];

for k=1:N-1
    [R_true(:,:,k+1) W_true(:,k+1)]=lgviSO(R_true(:,:,k),W_true(:,k));
end

%% Sampling attitude and angular velocity measurements
disp('Generating the measurements...');

rng(1);
H=diag([1.8 1.6 2.4]);
W_err_sigma=h*H*H';
%W_err_sigma=h*diag([0.5^2 0.8^2 1^2]);
W_err = (randn(N,3)*chol(W_err_sigma))';
W_mea = W_true + W_err;

Fz=diag([40 50 35]);
[R_mea_err accept_ratio]=pdf_MF_sampling(Fz,N);
for k=1:N
    R_mea(:,:,k)=R_true(:,:,k)*R_mea_err(:,:,k);
    rot_mea_err(k)=180/pi*norm(vee(logm(R_mea_err(:,:,k))));
end


%% Estimation routine
disp('Estimating...');
F=zeros(3,3,N); U=zeros(3,3,N); S=zeros(3,3,N); V=zeros(3,3,N); RM=zeros(3,3,N);

F(:,:,1)=F0;

[U0 S0 V0]=psvd(F0);
U(:,:,1)=U0;
S(:,:,1)=S0;
V(:,:,1)=V0;
RM(:,:,1)=U(:,:,1)*V(:,:,1)';

k_R_mea=[];

t_start=tic;
for k=1:N-1    
    % Prediction / Propagation
    if EST_METHOD == 0
        [F(:,:,k+1) U(:,:,k+1) S(:,:,k+1) V(:,:,k+1)]=FirstOrderPropagation(F(:,:,k),W_mea(:,k),W_err_sigma,h);
    else
        [F(:,:,k+1) U(:,:,k+1) S(:,:,k+1) V(:,:,k+1)]=UnscentedPropagation(F(:,:,k),W_mea(:,k),W_err_sigma,h);
    end
    RM(:,:,k+1)=U(:,:,k+1)*V(:,:,k+1)';        
    
    if rem(k,5)==0         
        % Measurement update / Correction
        k_R_mea=[k_R_mea; k+1];
        [F(:,:,k+1) U(:,:,k+1) S(:,:,k+1) V(:,:,k+1) RM(:,:,k+1)]=MeasurementUpdate(F(:,:,k+1),R_mea(:,:,k+1),Fz);
    end    
    
    if rem(k,5)==0
        disp(['Completed ' num2str(k/(N-1)*100,3) '%']);
        disp([F(:,:,k+1), S(:,:,k+1)]);
        disp([180/pi*norm((logmso3(R_true(:,:,k)'*RM(:,:,k))))]);
        disp(' ');        
    end    
end
t_elapsed = toc(t_start);

disp('Post processing...');
for k=1:N
    errRf(k)=norm(R_true(:,:,k)-RM(:,:,k));
    rot_est_err(k)=180/pi*norm((logmso3(R_true(:,:,k)'*RM(:,:,k))));
    s(:,k)=diag(S(:,:,k));
end

figure;
plot(t,rot_est_err,'b',t(k_R_mea),rot_est_err(k_R_mea),'bo');
xlabel('$t$','interpreter','latex');
ylabel('Est. error');
figure;
for ii=1:3
    subplot(3,1,ii);
    plot(t,1./s(ii,:),'b',t(k_R_mea),1./s(ii,k_R_mea),'b.');
end
subplot(3,1,1);ylabel('$\frac{1}{s_2+s_3}$','interpreter','latex');
subplot(3,1,2);ylabel('$\frac{1}{s_3+s_1}$','interpreter','latex');
subplot(3,1,3);ylabel('$\frac{1}{s_1+s_2}$','interpreter','latex');
xlabel('$t$','interpreter','latex');

%%
save(filename);
evalin('base',['load ' filename]);
disp(['Saved to "' filename '.mat", and loaded to workspace']);

end

function [F U S V RM]=MeasurementUpdate(F,Z,Fz)
F=F+Z*Fz';

[U S V]=psvd(F);
RM=U*V';
end

function [Fkp Ukp Skp Vkp]=FirstOrderPropagation(Fk,Wk,Gk,h)
[Uk Sk Vk]=psvd(Fk);
EQk=pdf_MF_moment(diag(Sk));
ERk=Uk*diag(EQk)*Vk';

ERkp=ERk*(eye(3)+h/2*(-trace(Gk)*eye(3)+Gk))*expm(h*hat(Wk));

[Ukp Dkp Vkp]=psvd(ERkp);

skp=pdf_MF_M2S(diag(Dkp));
Skp=diag(skp);

Fkp=Ukp*Skp*Vkp';
end

function [Fkp Ukp Skp Vkp]=UnscentedPropagation(Fk,Wk,Gk,h)

[R W0 W]=pdf_MF_unscented_transform(Fk);

R_bar_kp=W0*R(:,:,1)*expm(hat(h*Wk));
for i=1:3
    R_bar_kp=R_bar_kp+W(i)*R(:,:,2*i)*expm(hat(h*Wk))+W(i)*R(:,:,2*i+1)*expm(hat(h*Wk));
end
R_bar_kp=R_bar_kp*(eye(3)+h/2*(-trace(Gk)*eye(3)+Gk));

[Ukp Dkp Vkp]=psvd(R_bar_kp);

skp=pdf_MF_M2S(diag(Dkp));
Skp=diag(skp);

Fkp=Ukp*Skp*Vkp';
end



function [Rkp Wkp]=lgviSO(Rk,Wk)
global m g rho J e1 e2 e3 h

Mgk=m*g*hat(rho)*Rk'*e3;
fi=h*Wk+h^2/2*inv(J)*Mgk;
f=fi+1;
gk=h*J*Wk+h^2/2*Mgk;
while norm(fi-f) > 1e-15
    f=fi;
    GG=gk+(hat(gk)-2*J)*f+(gk'*f)*f;
    nabG=(hat(gk)-2*J)+gk'*f*eye(3)+f*gk';
    fi=f-inv(nabG)*GG;
end
Fk=(eye(3)+hat(fi))*inv(eye(3)-hat(fi));
Rkp=Rk*Fk;
Mgkp=m*g*hat(rho)*Rkp'*e3;
Wkp=J\(Fk'*(J*Wk+h/2*Mgk)+h/2*Mgkp);

end
