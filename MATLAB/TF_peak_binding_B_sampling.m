function B_new = TF_peak_binding_B_sampling(A, B, T, B_state, B_mean, B_var, sigma_A_noise)
% B_new=B;
M=size(B,2);
S=size(A,2);
F=size(A,1);
TF_index=randperm(M);

for i=1:M
    m=TF_index(i);
    temp_var=T(m,:)*T(m,:)'*B_var(m)/S+sigma_A_noise;%1 by 1 number
    mean_B=((A-B*T+B(:,m)*T(m,:))*T(m,:)'*B_var(m)/S+B_mean(:,m)*sigma_A_noise)/temp_var;% F by 1 vector
    vairance_B=B_var(m)*sigma_A_noise/temp_var;%1 by 1 number
    B(:,m)=(randn(F,1)*sqrt(vairance_B)+mean_B).*B_state(:,m);% F by 1 vector, conditional on B_state(:,t)
end
B_new=B;