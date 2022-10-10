function  [Gene_names, Gene_log2FC, Gene_pvalue, Gene_robust_pvalue]=RNA_psuedo_bulk_differential(Gene_diff_symbol, RNA, Case_sample_idx, Control_sample_idx)
RNA_samples = RNA.Properties.VariableNames(2:end);
RNA_count = RNA{:,2:end};
RNA_symbols=RNA.Var1;
RNA_flag=RNA_count;
RNA_flag(RNA_flag==1)=0;
RNA_flag(RNA_flag>0)=1;
total_RNA_reads=5e6;
RNA_raw_count_sum=sum(RNA_count)+1;

for s=1:length(RNA_samples)
    RNA_count(:,s)=RNA_count(:,s)/RNA_raw_count_sum(s)*total_RNA_reads;
end


RNA_count=RNA_count(:, [Case_sample_idx, Control_sample_idx]);
RNA_flag=RNA_flag(:, [Case_sample_idx, Control_sample_idx]);
Sample_label=[ones(1,length(Case_sample_idx)), -1*ones(1,length(Control_sample_idx))];


[Gene_names, gene_idx]=intersect(RNA_symbols, intersect(RNA_symbols(sum(RNA_flag,2)>=2), Gene_diff_symbol));
Gene_log2count=log2(RNA_count(gene_idx,:)+1);
Gene_pvalue=-1*ones(length(gene_idx), 1);
Gene_robust_pvalue=-1*ones(length(gene_idx), 1);
Gene_log2FC=zeros(length(gene_idx), 1);

for i=1:length(Gene_names)
    Gene_log2FC(i)=mean(Gene_log2count(i,Sample_label==1))-mean(Gene_log2count(i,Sample_label==-1));
    
    mdl = fitlm(Sample_label', Gene_log2count(i,:)');
    Gene_pvalue(i) = table2array(mdl.Coefficients(2,4));
    
    [B,STATS] = robustfit(Sample_label', Gene_log2count(i,:)');
    Gene_robust_pvalue(i) = STATS.p(2);
end

