function [Common_samples, Candidate_TFs, Candidate_TF_log2Count,...
    Candidate_Peaks, Candidate_Peak_log2Count,...
    Candidate_Genes, Candidate_Gene_TSS, Candidate_Gene_log2Count,...
    Candidate_TF_Peak_Binding, Candidate_Peak_Gene_looping]=  ...
    Initial_integration(Diff_Peaks, Peaks_bulk, ATAC_bulk, ATAC_Case_samples, ATAC_Control_samples, ...
                        Diff_Genes, RNA_bulk, RNA_Case_samples, RNA_Control_samples,...
                        TFs, TF_peaks, TF_Peak_Binding)

                    
                    
                    
%% select common samples between ATAC and RNA
ATAC_sample_names = ATAC_bulk.Properties.VariableNames(2:end);
ATAC_count = ATAC_bulk{:,2:end};
[ATAC_sample_names, sample_index]=intersect(ATAC_sample_names, union(ATAC_Case_samples, ATAC_Control_samples));
ATAC_count=ATAC_count(:,sample_index);


RNA_symbols = RNA_bulk.Var1;
RNA_sample_names = RNA_bulk.Properties.VariableNames(2:end);
RNA_count = RNA_bulk{:,2:end};
[RNA_sample_names, sample_index]=intersect(RNA_sample_names, union(RNA_Case_samples, RNA_Control_samples));
RNA_count=RNA_count(:,sample_index);



[Common_samples, bidex, cidex] = intersect(ATAC_sample_names, RNA_sample_names);
ATAC_count = ATAC_count(:,bidex);
RNA_count = RNA_count(:,cidex);



%% bulk data normalization, bulk data will be used in vairiable prior calculation
%we take out the mean of all samples (assumed to be baseline activity)
%we don't take z-score as peaks or genes with high varaince across samples between two groups
%always mean more signal than low variance ones

total_ATAC_reads=5e6;
ATAC_raw_count_sum=sum(ATAC_count)+1;
for s=1:length(Common_samples)
    ATAC_count(:,s)=ATAC_count(:,s)/ATAC_raw_count_sum(s)*total_ATAC_reads;
end

[Candidate_Peaks, bidex, cidex]=intersect(Diff_Peaks, Peaks_bulk, 'stable', 'rows');
Candidate_Peak_log2Count=log2(ATAC_count(cidex,:)+1);
for p=1:size(Candidate_Peaks, 1)
    Candidate_Peak_log2Count(p,:)=Candidate_Peak_log2Count(p,:)-mean(Candidate_Peak_log2Count(p,:));
end

[Candidate_Peaks, bidex, cidex]=intersect(Candidate_Peaks, TF_peaks, 'stable', 'rows');
Candidate_Peak_log2Count=Candidate_Peak_log2Count(bidex,:);
Candidate_TF_Peak_Binding=TF_Peak_Binding(cidex,:);


total_RNA_reads=1e6;
RNA_raw_count_sum=sum(RNA_count)+1;
for s=1:length(Common_samples)
    RNA_count(:,s)=RNA_count(:,s)/RNA_raw_count_sum(s)*total_RNA_reads;
end

[Candidate_Genes, bidex, cidex]=intersect(Diff_Genes, RNA_symbols, 'stable');

Candidate_Gene_log2Count=log2(RNA_count(cidex,:)+1);
for g=1:length(Candidate_Genes)
    Candidate_Gene_log2Count(g,:)=Candidate_Gene_log2Count(g,:)-mean(Candidate_Gene_log2Count(g,:));
end




%% pair peaks to genes within 500k bps, refered to RefSeq hg38 file, this can be replaced using HiC data

[hg38.chr, hg38.strand, hg38.start, hg38.end, hg38.gene_name]=textread('hg38_Refseq', '%d %s %d %d %s', 'headerlines', 1);
[Candidate_Genes, bidex]=intersect(Candidate_Genes, hg38.gene_name);
Candidate_Gene_log2Count=Candidate_Gene_log2Count(bidex,:);
Candidate_Gene_TSS=zeros(length(Candidate_Genes), 2);% chr, TSS

for g=1:length(Candidate_Genes)
    index=find(strcmp(hg38.gene_name, Candidate_Genes{g})>0);
    if ~isempty(index)
        Candidate_Gene_TSS(g,1)=hg38.chr(index(1));
        if strcmp(hg38.strand(index(1)), '+')>0
            Candidate_Gene_TSS(g,2)=min(hg38.start(index));
        else
            Candidate_Gene_TSS(g,2)=max(hg38.end(index));
        end
    end
end


Candidate_Peak_Gene_looping=sparse(length(Candidate_Peaks), length(Candidate_Genes));
for g=1:length(Candidate_Genes)
    index=find(Candidate_Peaks(:,1)-Candidate_Gene_TSS(g,1)==0 & abs((Candidate_Peaks(:,2)+Candidate_Peaks(:,3))/2-Candidate_Gene_TSS(g,2))<5e5);
    if ~isempty(index)
        Candidate_Peak_Gene_looping(index,g)=1;
    end
end



%% select actively expressed TFs
[Candidate_TFs, bidex, cidex]=intersect(TFs, RNA_symbols, 'stable');
Candidate_TF_Peak_Binding=Candidate_TF_Peak_Binding(:,bidex);
TF_Peak_Binding=TF_Peak_Binding(:,bidex);
Candidate_TF_log2Count=log2(RNA_count(cidex,:)+1);

active_TF_index=find(sum(Candidate_TF_log2Count,2)>0);
Candidate_TFs=Candidate_TFs(active_TF_index);
TF_Peak_Binding=TF_Peak_Binding(:,active_TF_index);
Candidate_TF_Peak_Binding=Candidate_TF_Peak_Binding(:,active_TF_index);
Candidate_TF_log2Count=Candidate_TF_log2Count(active_TF_index,:);



%% After initialization, TF, peaks and genes and all their priors 
Peak_index=find(sum(Candidate_Peak_Gene_looping,2)>0);
Candidate_Peaks=Candidate_Peaks(Peak_index,:);
Candidate_Peak_log2Count=Candidate_Peak_log2Count(Peak_index,:);

Gene_index=find(sum(Candidate_Peak_Gene_looping)>0);
Candidate_Genes=Candidate_Genes(Gene_index);
Candidate_Gene_TSS=Candidate_Gene_TSS(Gene_index,:);
Candidate_Gene_log2Count=Candidate_Gene_log2Count(Gene_index,:);

Candidate_Peak_Gene_looping=Candidate_Peak_Gene_looping(Peak_index,Gene_index);


TF_index=find(sum(Candidate_TF_Peak_Binding(Peak_index,:))/length(Peak_index)>0.1 &...
    sum(full(Candidate_TF_Peak_Binding(Peak_index,:)))/length(Peak_index)./(full(sum(TF_Peak_Binding))/length(TF_peaks))>1.2);

Candidate_TF_Peak_Binding=Candidate_TF_Peak_Binding(Peak_index,TF_index);
Candidate_TFs=Candidate_TFs(TF_index);
Candidate_TF_log2Count=Candidate_TF_log2Count(TF_index,:);



