function [T_A_prior, T_R_prior, T_sample_mean, T_sample_var, B_prior, B_mean, B_var, B_prob, L_prior, L_mean, L_var, L_prob]=MAGICAL_initialization(Candidate_TF_log2Count,...
        Candidate_Peak_log2Count, Candidate_Gene_log2Count, Candidate_TF_Peak_Binding, Candidate_Peak_Gene_looping, M, S, ATAC_cell_vector, RNA_cell_vector)

% TF activity prior: TF RNA expression

T_sample_mean=Candidate_TF_log2Count;
T_sample_var=var(Candidate_TF_log2Count,0,2);

T_A_prior=zeros(M, length(ATAC_cell_vector));
T_R_prior=zeros(M, length(RNA_cell_vector));
for m=1:M
    for s=1:S
        A_index=find(ATAC_cell_vector==s);
        T_A_prior(m,A_index)=randn(1,length(A_index))*sqrt(T_sample_var(m))+T_sample_mean(m,s);
        R_index=find(RNA_cell_vector==s);
        T_R_prior(m,R_index)=randn(1,length(R_index))*sqrt(T_sample_var(m))+T_sample_mean(m,s);
    end
end



% TF-peak binding prior: regression weight between TF activity (expression) and peak ATAC activity

B_prior=full(Candidate_TF_Peak_Binding);%mean per binding variable
B_prob=full(Candidate_TF_Peak_Binding);
B_var=ones(M,1);
[xx,yy]=find(Candidate_TF_Peak_Binding>0);
for i=1:length(xx)
    mdl = fitlm(Candidate_TF_log2Count(yy(i),:)', Candidate_Peak_log2Count(xx(i),:)');
    B_prior(xx(i), yy(i)) = table2array(mdl.Coefficients(2,1)); 
    B_prob(xx(i), yy(i)) = 1-table2array(mdl.Coefficients(2,4));
end
B_mean=B_prior;
for m=1:M
    B_var(m)=1*var(B_prior(Candidate_TF_Peak_Binding(:,m)>0,m));%variance of binding events weights of TF t
end


% Peak-Gene looping prior

L_prior=full(Candidate_Peak_Gene_looping);
L_prob=full(Candidate_Peak_Gene_looping);
[xx,yy]=find(Candidate_Peak_Gene_looping>0);
for i=1:length(xx)
    mdl = fitlm(Candidate_Peak_log2Count(xx(i),:)', Candidate_Gene_log2Count(yy(i),:)');
    ATAC_RNA_weight_1 = table2array(mdl.Coefficients(2,1));
    ATAC_RNA_probability_1 = 1-table2array(mdl.Coefficients(2,4));
    L_prior(xx(i),yy(i))=ATAC_RNA_weight_1;
    L_prob(xx(i),yy(i))=ATAC_RNA_probability_1;
end
L_mean=L_prior;
L_var=var(L_prior(L_prior~=0));
