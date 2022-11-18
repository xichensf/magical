function     B = TF_peak_binding_B_sampling(scATAC_read_count_matrix, ATAC_cell_vector, A_sample, B, T_A, B_state, B_mean, B_var, sigma_A_noise, M, S, P, G)


T_sample=zeros(M,S);
for s=1:S
    T_sample(:,s)=mean(T_A(:,ATAC_cell_vector==s), 2);
end

TF_index=randperm(M);
for i=1:M
    m=TF_index(i);
    temp_var=T_sample(m,:)*T_sample(m,:)'*B_var(m)/S+sigma_A_noise;%1 by 1 number
    mean_B=((A_sample-B*T_sample+B(:,m)*T_sample(m,:))*T_sample(m,:)'*B_var(m)/S+B_mean(:,m)*sigma_A_noise)/temp_var;% P by 1 vector
    vairance_B=B_var(m)*sigma_A_noise/temp_var;%1 by 1 number
    bb=randn(P,1);
    bb(bb>3)=3;
    bb(bb<-3)=-3;
    B(:,m)=(bb*sqrt(vairance_B)+mean_B).*B_state(:,m);% P by 1 vector, conditional on B_state(:,t)
end