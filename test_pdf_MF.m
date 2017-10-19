clear all;
close all;

F=rand(3,3);
[U S V]=psvd(F);
s=diag(S);
F=U*diag(s)*V';

c=pdf_MF_normal(s);
c_bar=pdf_MF_normal(s,1);

disp(['check:c_bar    ' num2str(norm(abs(c-c_bar*exp(sum(s)))/c))]);

[dc ddc]=pdf_MF_normal_deriv(s,1,0);
[dc_bar ddc_bar]=pdf_MF_normal_deriv(s,1,1);
dc_new=exp(sum(s))*(dc_bar+c_bar*ones(3,1));
disp(['check:dc_bar   ' num2str(norm(dc-dc_new)/norm(dc))]);

ddc_new=zeros(3,3);
for i=1:3
    for j=1:3
        ddc_new(i,j)=c_bar+dc_bar(i)+dc_bar(j)+ddc_bar(i,j);
    end
end
ddc_new=ddc_new*exp(sum(s));

disp(['check:ddc_bar  ' num2str(norm(ddc-ddc_new)/norm(ddc))]);

[M1, M2]=pdf_MF_moment(s,1);
disp(['check:M1       ' num2str(norm(M1-dc/c))]);
disp(['check:M2       ' num2str(norm(M2-ddc/c))]);

[s_new f_final niter]=pdf_MF_M2S(M1);
disp(['check:M2S      ' num2str(norm(s-s_new))]);

[R W0 W theta sigma_min sigma]=pdf_MF_unscented_transform(F);
M1_sigma=W0*R(:,:,1);
for i=1:3
    M1_sigma=M1_sigma+W(i)*R(:,:,2*i)+W(i)*R(:,:,2*i+1);
end

disp(['check:UT       ' num2str(norm(M1_sigma-U*diag(M1)*V'))]);

F_new=pdf_MF_inv_unscented_transform(R,W0,W);
disp(['check:invUT    ' num2str(norm(F-F_new))]);




