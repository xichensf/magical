function [Candidate_TFs, Candidate_TF_log2Count,...
    Candidate_peaks, Candidate_Peak_log2Count,...
    Candidate_genes, Candidate_Gene_log2Count,...
    Candidate_TF_Peak_Binding, Candidate_Peak_Gene_looping]=...
    Candidate_circuits_construction_with_TAD(Common_samples, Candidate_genes, Candidate_peaks,...
                                  scRNA_genes, scRNA_cells, scRNA_read_count_matrix, ...
                                  scATAC_peaks, scATAC_cells, scATAC_read_count_matrix,...
                                  Motifs, TF_peak_binding_matrix,...
                                  Refseq, TAD)

fprintf('Candidate regulatory circuits constrcution ...\n\n')

% TF-peak binding
[aa, bidex, cidex]=intersect([scATAC_peaks.chr_num, scATAC_peaks.point1, scATAC_peaks.point2], [Candidate_peaks.chr_num, Candidate_peaks.point1, Candidate_peaks.point2], 'rows', 'stable');
Candidate_peaks.peak_index=scATAC_peaks.peak_index(bidex);
Candidate_TF_Peak_Binding=TF_peak_binding_matrix(bidex,:);

Candidate_peaks.chr=Candidate_peaks.chr(cidex);
Candidate_peaks.chr_num=Candidate_peaks.chr_num(cidex);
Candidate_peaks.point1=Candidate_peaks.point1(cidex);
Candidate_peaks.point2=Candidate_peaks.point2(cidex);


TF_num=full(sum(Candidate_TF_Peak_Binding));
TF_pct=full(sum(Candidate_TF_Peak_Binding))/length(Candidate_peaks.peak_index);
TF_enrichment_FC=zeros(1,length(length(Motifs.name)));

for m=1:length(Motifs.name)
    TF_enrichment_FC(m)=(TF_num(m)/length(Candidate_peaks.peak_index))/(sum(TF_peak_binding_matrix(:,m))/length(scATAC_peaks.peak_index));
end
TF_index=find(TF_pct>0.05 & TF_num>30 & TF_enrichment_FC>1.2);

if ~isempty(TF_index)
    Candidate_TFs=Motifs.name(TF_index);
    Candidate_TF_Peak_Binding=Candidate_TF_Peak_Binding(:,TF_index);
else
    fprintf('Too few peaks with TF binding sites.\nMAGICAL not applicable to this cell type!\n\n\n\n')
    %         continue;
end

peak_index=find(sum(Candidate_TF_Peak_Binding, 2)>0);
Candidate_peaks.peak_index=Candidate_peaks.peak_index(peak_index);
Candidate_peaks.chr=Candidate_peaks.chr(peak_index);
Candidate_peaks.chr_num=Candidate_peaks.chr_num(peak_index);
Candidate_peaks.point1=Candidate_peaks.point1(peak_index);
Candidate_peaks.point2=Candidate_peaks.point2(peak_index);
Candidate_TF_Peak_Binding=Candidate_TF_Peak_Binding(peak_index,:);


% Peak-gene looping
Candidate_genes.gene_symbols=intersect(Candidate_genes.gene_symbols, Refseq.gene_name);
Candidate_genes.gene_TSS=zeros(length(Candidate_genes.gene_symbols), 2);% chr, TSS

for g=1:length(Candidate_genes.gene_symbols)
    index=find(strcmp(Refseq.gene_name, Candidate_genes.gene_symbols{g})>0);
    if ~isempty(index)
        Candidate_genes.gene_TSS(g,1)=Refseq.chr(index(1));
        if strcmp(Refseq.strand(index(1)), '+')>0
            Candidate_genes.gene_TSS(g,2)=min(Refseq.start(index));
        else
            Candidate_genes.gene_TSS(g,2)=max(Refseq.end(index));
        end
    end
end


Candidate_Peak_Gene_looping=sparse(length(Candidate_peaks.peak_index), length(Candidate_genes.gene_symbols));
for t=1:length(TAD.chr_num)
    if TAD.right_boundary(t)-TAD.left_boundary(t)<2e6
        peak_index=find(Candidate_peaks.chr_num==TAD.chr_num(t) & Candidate_peaks.point1>TAD.left_boundary(t) & Candidate_peaks.point2<TAD.right_boundary(t));
        gene_index=find(Candidate_genes.gene_TSS(:,1)==TAD.chr_num(t) & Candidate_genes.gene_TSS(:,2)>TAD.left_boundary(t) & Candidate_genes.gene_TSS(:,2)<TAD.right_boundary(t));
        if ~isempty(peak_index) && ~isempty(gene_index)
            Candidate_Peak_Gene_looping(peak_index,gene_index)=1;
        end
    end
end


peak_index=find(sum(Candidate_Peak_Gene_looping,2)>0);
gene_index=find(sum(Candidate_Peak_Gene_looping)>0);
Candidate_Peak_Gene_looping=Candidate_Peak_Gene_looping(peak_index,gene_index);

Candidate_peaks.peak_index=Candidate_peaks.peak_index(peak_index);
Candidate_peaks.chr=Candidate_peaks.chr(peak_index);
Candidate_peaks.chr_num=Candidate_peaks.chr_num(peak_index);
Candidate_peaks.point1=Candidate_peaks.point1(peak_index);
Candidate_peaks.point2=Candidate_peaks.point2(peak_index);

Candidate_genes.gene_symbols=Candidate_genes.gene_symbols(gene_index);
Candidate_genes.gene_TSS=Candidate_genes.gene_TSS(gene_index,:);% chr, TSS

TF_index=find(sum(Candidate_TF_Peak_Binding(peak_index,:))>0);
Candidate_TFs=Candidate_TFs(TF_index);
Candidate_TF_Peak_Binding=Candidate_TF_Peak_Binding(peak_index, TF_index);



%filter the input paired data for candidate peaks and candidate genes

for s=1:length(Common_samples)
    ATAC_count(:,s)=full(sum(scATAC_read_count_matrix(:,strcmp(scATAC_cells.subject_ID, Common_samples{s})), 2));
    RNA_count(:,s)=full(sum(scRNA_read_count_matrix(:,strcmp(scRNA_cells.subject_ID, Common_samples{s})), 2));
end


total_ATAC_reads=5e6;
ATAC_raw_count_sum=sum(ATAC_count)+1;
for s=1:length(Common_samples)
    ATAC_count(:,s)=ATAC_count(:,s)/ATAC_raw_count_sum(s)*total_ATAC_reads;
end
ATAC_log2=log2(ATAC_count+1);


total_RNA_reads=5e6;
RNA_raw_count_sum=sum(RNA_count)+1;
for s=1:length(Common_samples)
    RNA_count(:,s)=RNA_count(:,s)/RNA_raw_count_sum(s)*total_RNA_reads;
end
RNA_log2=log2(RNA_count+1);


% select actively accessbile peaks
[aa, bidex, cidex]=intersect(Candidate_peaks.peak_index, scATAC_peaks.peak_index, 'stable');
Candidate_peaks.peak_index=Candidate_peaks.peak_index(bidex);
Candidate_peaks.chr=Candidate_peaks.chr(bidex);
Candidate_peaks.chr_num=Candidate_peaks.chr_num(bidex);
Candidate_peaks.point1=Candidate_peaks.point1(bidex);
Candidate_peaks.point2=Candidate_peaks.point2(bidex);
Candidate_TF_Peak_Binding=Candidate_TF_Peak_Binding(bidex,:);
Candidate_Peak_Gene_looping=Candidate_Peak_Gene_looping(bidex,:);
Candidate_Peak_log2Count=ATAC_log2(cidex,:);

peak_index=find(sum(Candidate_Peak_log2Count,2)>0);
Candidate_peaks.peak_index=Candidate_peaks.peak_index(peak_index);
Candidate_peaks.chr=Candidate_peaks.chr(peak_index);
Candidate_peaks.chr_num=Candidate_peaks.chr_num(peak_index);
Candidate_peaks.point1=Candidate_peaks.point1(peak_index);
Candidate_peaks.point2=Candidate_peaks.point2(peak_index);
Candidate_TF_Peak_Binding=Candidate_TF_Peak_Binding(peak_index,:);
Candidate_Peak_Gene_looping=Candidate_Peak_Gene_looping(peak_index,:);
Candidate_Peak_log2Count=Candidate_Peak_log2Count(peak_index,:);

for p=1:length(Candidate_peaks.peak_index)
    Candidate_Peak_log2Count(p,:)=Candidate_Peak_log2Count(p,:)-mean(Candidate_Peak_log2Count(p,:));
end


% select actively expressed genes
[aa, bidex, cidex]=intersect(Candidate_genes.gene_symbols, scRNA_genes.gene_symbols, 'stable');
Candidate_genes.gene_symbols=Candidate_genes.gene_symbols(bidex);
Candidate_genes.gene_TSS=Candidate_genes.gene_TSS(bidex,:);
Candidate_Peak_Gene_looping=Candidate_Peak_Gene_looping(:,bidex);
Candidate_Gene_log2Count=RNA_log2(cidex,:);

gene_index=find(sum(Candidate_Gene_log2Count,2)>0);
Candidate_genes.gene_symbols=Candidate_genes.gene_symbols(gene_index);
Candidate_genes.gene_TSS=Candidate_genes.gene_TSS(gene_index,:);
Candidate_Peak_Gene_looping=Candidate_Peak_Gene_looping(:,gene_index);
Candidate_Gene_log2Count=Candidate_Gene_log2Count(gene_index,:);

for g=1:length(Candidate_genes.gene_symbols)
    Candidate_Gene_log2Count(g,:)=Candidate_Gene_log2Count(g,:)-mean(Candidate_Gene_log2Count(g,:));
end


% select actively expressed TFs
[Candidate_TFs, bidex, cidex]=intersect(Candidate_TFs, scRNA_genes.gene_symbols, 'stable');
Candidate_TF_Peak_Binding=Candidate_TF_Peak_Binding(:,bidex);
Candidate_TF_log2Count=RNA_log2(cidex,:);

active_TF_index=find(sum(Candidate_TF_log2Count,2)>0);
Candidate_TFs=Candidate_TFs(active_TF_index);
Candidate_TF_Peak_Binding=Candidate_TF_Peak_Binding(:,active_TF_index);
Candidate_TF_log2Count=Candidate_TF_log2Count(active_TF_index,:);


%final filter and candidate circuits selection
peak_index=find(sum(Candidate_Peak_Gene_looping,2)>0 & sum(Candidate_TF_Peak_Binding,2)>0);
gene_index=find(sum(Candidate_Peak_Gene_looping(peak_index,:))>0);
TF_index=find(sum(Candidate_TF_Peak_Binding(peak_index,:))>0);

Candidate_TF_Peak_Binding=Candidate_TF_Peak_Binding(peak_index,TF_index);
Candidate_Peak_Gene_looping=Candidate_Peak_Gene_looping(peak_index,gene_index);


Candidate_TFs=Candidate_TFs(TF_index);
Candidate_TF_log2Count=Candidate_TF_log2Count(TF_index,:);


Candidate_genes.gene_symbols=Candidate_genes.gene_symbols(gene_index);
Candidate_genes.gene_TSS=Candidate_genes.gene_TSS(gene_index,:);
Candidate_Gene_log2Count=Candidate_Gene_log2Count(gene_index,:);


Candidate_peaks.peak_index=Candidate_peaks.peak_index(peak_index);
Candidate_peaks.chr=Candidate_peaks.chr(peak_index);
Candidate_peaks.chr_num=Candidate_peaks.chr_num(peak_index);
Candidate_peaks.point1=Candidate_peaks.point1(peak_index);
Candidate_peaks.point2=Candidate_peaks.point2(peak_index);
Candidate_Peak_log2Count=Candidate_Peak_log2Count(peak_index,:);
