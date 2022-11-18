function  L = Peak_gene_looping_L_samping(scRNA_read_count_matrix, RNA_cell_vector, R_sample, L, B, T_R, L_state, L_mean, L_var, sigma_R_noise, M, S, P, G)

T_sample=zeros(M,S);
for s=1:S
    T_sample(:,s)=mean(T_R(:,RNA_cell_vector==s), 2);
end
  

A_estimate=B*T_sample;

P_index=randperm(P);
for i=1:P
    p=P_index(i);
    temp_var=A_estimate(p,:)*A_estimate(p,:)'*L_var/S+sigma_R_noise;%1 by 1 number
    mean_L=(A_estimate(p,:)*(R_sample-L'*A_estimate+L(p,:)'*A_estimate(p,:))'*L_var/S+L_mean(p,:)*sigma_R_noise)/temp_var;%1 by G vector
    vairance_L=L_var*sigma_R_noise/temp_var;%1 by 1 number
    ll=randn(1,G);
    ll(ll>3)=3;
    ll(ll<-3)=-3;
    L(p,:)=(ll*sqrt(vairance_L)+mean_L).*L_state(p,:);%1 by G vector
end
