function T_new = TF_activity_T_sampling(A, B, T, T_mean, T_var, sigma_A_noise)

M=size(T,1);
S=size(A,2);
F=size(A,1);
TF_index=randperm(M);

for i=1:M
    m=TF_index(i);
    temp_var=(B(:,m)'*B(:,m)*T_var(m)/F+sigma_A_noise);% 1 by 1 number
    mean_T=(B(:,m)'*(A-B*T+B(:,m)*T(m,:))/F+T_mean(m,:)*sigma_A_noise)/temp_var;%1 by S vector
    variance_T=T_var(m)*sigma_A_noise/temp_var;%1 by 1 number
    aa=randn(1,S);
    aa(aa>3)=3;
    aa(aa<-3)=-3;
    T(m,:)=aa*sqrt(variance_T)+mean_T;%1 by S vector
end
T_new=T;