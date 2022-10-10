function [Candidate_Peaks, ATAC_log2FC, ATAC_pvalue, ATAC_robust_pvalue]=ATAC_psuedo_bulk_differential(Candidate_Peaks, Peaks_bulk, ATAC, Case_sample_idx, Control_sample_idx)
ATAC_samples = ATAC.Properties.VariableNames(2:end);
ATAC_count = ATAC{:,2:end};
ATAC_flag = ATAC_count;
ATAC_flag(ATAC_flag==1) = 0;
ATAC_flag(ATAC_flag>0) = 1;
total_ATAC_reads=5e6;
ATAC_raw_count_sum=sum(ATAC_count)+1;
for s=1:length(ATAC_samples)
    ATAC_count(:,s)=ATAC_count(:,s)/ATAC_raw_count_sum(s)*total_ATAC_reads;
end

ATAC_count=ATAC_count(:, [Case_sample_idx, Control_sample_idx]);
ATAC_flag=ATAC_flag(:, [Case_sample_idx, Control_sample_idx]);
Sample_label=[ones(1,length(Case_sample_idx)), -1*ones(1,length(Control_sample_idx))];

%active peaks in at least two samples
[Candidate_Peaks, peak_idx]=intersect(Peaks_bulk, intersect(Peaks_bulk(sum(ATAC_flag,2)>=2,:), Candidate_Peaks, 'rows'), 'rows', 'stable');
ATAC_log2count=log2(ATAC_count(peak_idx,:)+1);
ATAC_log2FC=zeros(length(peak_idx), 1);
ATAC_pvalue=-1*ones(length(peak_idx), 1);
ATAC_robust_pvalue=-1*ones(length(peak_idx), 1);

for i=1:size(Candidate_Peaks, 1)
    ATAC_log2FC(i)=mean(ATAC_log2count(i,Sample_label==1))-mean(ATAC_log2count(i,Sample_label==-1));
    
    mdl = fitlm(Sample_label', ATAC_log2count(i,:)');
    ATAC_pvalue(i) = table2array(mdl.Coefficients(2,4));
    
    [B,STATS] = robustfit(Sample_label', ATAC_log2count(i,:)');
    ATAC_robust_pvalue(i) = STATS.p(2);
end
