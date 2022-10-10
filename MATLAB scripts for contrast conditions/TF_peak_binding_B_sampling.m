function B_new = TF_peak_binding_B_sampling(A, B, P, B_state, B_mean, B_var, sigma_A_noise)
% B_new=B;
T=size(B,2);
S=size(A,2);
F=size(A,1);
T_index=randperm(T);

for i=1:T
    t=T_index(i);
    temp_var=P(t,:)*P(t,:)'*B_var(t)/S+sigma_A_noise;%1 by 1 number
    mean_B=((A-B*P+B(:,t)*P(t,:))*P(t,:)'*B_var(t)/S+B_mean(:,t)*sigma_A_noise)/temp_var;% F by 1 vector
    vairance_B=B_var(t)*sigma_A_noise/temp_var;%1 by 1 number
    B(:,t)=(randn(F,1)*sqrt(vairance_B)+mean_B).*B_state(:,t);% F by 1 vector, conditional on B_state(:,t)
end
B_new=B;