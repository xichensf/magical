function P_new = TF_activity_P_sampling(A, B, P, P_mean, P_var, sigma_A_noise)
% P_new=P;
T=size(P,1);
S=size(A,2);
F=size(A,1);
T_index=randperm(T);

for i=1:T
    t=T_index(i);
    temp_var=(B(:,t)'*B(:,t)*P_var(t)/F+sigma_A_noise);% 1 by 1 number
    mean_P=(B(:,t)'*(A-B*P+B(:,t)*P(t,:))/F+P_mean(t,:)*sigma_A_noise)/temp_var;%1 by S vector
    variance_P=P_var(t)*sigma_A_noise/temp_var;%1 by 1 number
    P(t,:)=randn(1,S)*sqrt(variance_P)+mean_P;%1 by S vector
end
P_new=P;