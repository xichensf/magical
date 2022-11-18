function [T_A, T_R, T_sample] = TF_activity_T_sampling(scATAC_read_count_matrix, ATAC_cell_vector, A_sample, scRNA_read_count_matrix, RNA_cell_vector, R_sample, B, T_A, T_R, T_mean, T_var, sigma_A_noise, M, S, P, G)

T_sample=zeros(M,S);

for s=1:S
    T_sample(:,s)=mean(T_A(:,ATAC_cell_vector==s),2);
end


TF_index=randperm(M);
for i=1:M
    m=TF_index(i);
    temp_var=(B(:,m)'*B(:,m)*T_var(m)/P+sigma_A_noise);% 1 by 1 number
    mean_T=(B(:,m)'*(A_sample-B*T_sample+B(:,m)*T_sample(m,:))/P+T_mean(m,:)*sigma_A_noise)/temp_var;%1 by S vector
    variance_T=T_var(m)*sigma_A_noise/temp_var;%1 by 1 number
    aa=randn(1,S);
    aa(aa>3)=3;
    aa(aa<-3)=-3;
    T_sample(m,:)=aa*sqrt(variance_T)+mean_T;%1 by S vector
    for s=1:S
        A_index=find(ATAC_cell_vector==s);
        T_A(m,A_index)=randn(1,length(A_index))*sqrt(variance_T)+T_sample(m,s);
        R_index=find(RNA_cell_vector==s);
        T_R(m,R_index)=randn(1,length(R_index))*sqrt(variance_T)+T_sample(m,s);        
    end
end
