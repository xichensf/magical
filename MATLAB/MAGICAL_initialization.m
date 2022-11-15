function [T_prior, T_mean, T_var, B_prior, B_mean, B_var, B_prob, L_prior, L_mean, L_var, L_prob]=MAGICAL_initialization(Candidate_TF_log2Count,...
        Candidate_Peak_log2Count, Candidate_Gene_log2Count, Candidate_TF_Peak_Binding, Candidate_Peak_Gene_looping, M)

% TF activity prior: TF RNA expression

T_prior=Candidate_TF_log2Count;%mean per TF activity per sample
T_mean=T_prior;
T_var=ones(M, 1);
for m=1:M
    T_var(m)=var(Candidate_TF_log2Count(m,:));%variance of protein activity of TF t
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
   
%     [B,STATS] = robustfit(Candidate_Peak_log2Count(xx(i),:)', Candidate_Gene_log2Count(yy(i),:)');
%     ATAC_RNA_weight_2 = B(2);
%     ATAC_RNA_probability_2 = 1-STATS.p(2);
%     L_prior(xx(i),yy(i))=mean([ATAC_RNA_weight_1, ATAC_RNA_weight_2]);
%     L_prob(xx(i),yy(i))=max(ATAC_RNA_probability_1, ATAC_RNA_probability_2);
end
L_mean=L_prior;
L_var=var(L_prior(L_prior~=0));
