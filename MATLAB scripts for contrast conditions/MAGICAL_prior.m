function [P_prior, P_mean, P_var, B_prior, B_mean, B_var, B_prob, L_prior, L_mean, L_var, L_prob]=MAGICAL_prior(Candidate_TF_log2Count,...
        Candidate_Peak_log2Count, Candidate_Gene_log2Count, Candidate_TF_Peak_Binding, Candidate_Peak_Gene_looping, S, T, F, G)


% TF activity prior: TF RNA expression

P_prior=Candidate_TF_log2Count;%mean per TF activity per sample
P_var=ones(T, 1);
for t=1:T
    P_prior(t,:)=Candidate_TF_log2Count(t,:)-mean(Candidate_TF_log2Count(t,:));
    P_var(t)=var(Candidate_TF_log2Count(t,:));%variance of protein activity of TF t
end
P_mean=P_prior;


% TF-peak binding prior: regression weight between TF activity (expression) and peak ATAC activity

B_prior=full(Candidate_TF_Peak_Binding);%mean per binding variable
B_prob=full(Candidate_TF_Peak_Binding);
B_var=ones(T,1);
[xx,yy]=find(Candidate_TF_Peak_Binding>0);
for i=1:length(xx)
    mdl = fitlm(Candidate_TF_log2Count(yy(i),:)', Candidate_Peak_log2Count(xx(i),:)');
    B_prior(xx(i), yy(i)) = table2array(mdl.Coefficients(2,1)); 
    B_prob(xx(i), yy(i)) = 1-table2array(mdl.Coefficients(2,4));
end
B_mean=B_prior;
for t=1:T
    B_var(t)=1*var(B_prior(Candidate_TF_Peak_Binding(:,t)>0,t));%variance of binding events weights of TF t
end


% Peak-Gene looping prior

L_prior=full(Candidate_Peak_Gene_looping);
L_prob=full(Candidate_Peak_Gene_looping);
[xx,yy]=find(Candidate_Peak_Gene_looping>0);
for i=1:length(xx)
    mdl = fitlm(Candidate_Peak_log2Count(xx(i),:)', Candidate_Gene_log2Count(yy(i),:)');
    ATAC_RNA_weight_1 = table2array(mdl.Coefficients(2,1));
    ATAC_RNA_probability_1 = 1-table2array(mdl.Coefficients(2,4));
    
    [B,STATS] = robustfit(Candidate_Peak_log2Count(xx(i),:)', Candidate_Gene_log2Count(yy(i),:)');
    ATAC_RNA_weight_2 = B(2);
    ATAC_RNA_probability_2 = 1-STATS.p(2);
    
    L_prior(xx(i),yy(i))=mean([ATAC_RNA_weight_1, ATAC_RNA_weight_2]);
    L_prob(xx(i),yy(i))=max(ATAC_RNA_probability_1, ATAC_RNA_probability_2);
end
L_mean=L_prior;
L_var=var(L_prior(L_prior~=0));