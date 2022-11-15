function MAGICAL(Candidate_gene_file_path, Candidate_peak_file_path,...
                scRNA_readcount_file_path, scRNA_gene_file_path, scRNA_cellmeta_file_path,...
                scATAC_readcount_file_path, scATAC_peak_file_path, scATAC_cellmeta_file_path,...
                Motif_mapping_file_path, Motif_name_file_path, TAD_flag, TAD_file_path, Ref_seq_file_path, ...
                Output_file_path, prob_threshold_TF_peak_binding, prob_threshold_peak_gene_looping, iteration_num)


fprintf('loading all input data ...\n\n')

%load candidate genes
Candidate_genes.gene_symbols=textread(Candidate_gene_file_path, '%s');


%load candidate peaks
[Candidate_peaks.chr, Candidate_peaks.point1, Candidate_peaks.point2]=textread(Candidate_peak_file_path, '%s %d %d');
Candidate_peaks.chr_num=zeros(length(Candidate_peaks.chr), 1);
for i=1:22
    Candidate_peaks.chr_num(strcmp(Candidate_peaks.chr, ['chr', num2str(i)])>0,1)=i;
    %we exclude X and Y chromosome by setting its chr index as 0
end


%load scRNAseq data
[scRNA_genes.gene_index, scRNA_genes.gene_symbols]=textread(scRNA_gene_file_path, '%d %s');

[scRNA_cells.cell_index, scRNA_cells.cell_barcode, scRNA_cells.cell_type, scRNA_cells.subject_ID, scRNA_cells.condition]=...
    textread(scRNA_cellmeta_file_path, '%d %s %s %s %s', 'delimiter', '\t');

[scRNA_read_count_table.gene_index, scRNA_read_count_table.cell_index, scRNA_read_count_table.readcount]=textread(scRNA_readcount_file_path, '%d %d %d');
scRNA_read_count_matrix=sparse(scRNA_read_count_table.gene_index,scRNA_read_count_table.cell_index,scRNA_read_count_table.readcount);
if size(scRNA_read_count_matrix,1)<length(scRNA_genes.gene_symbols)
    scRNA_read_count_matrix(end+1:length(scRNA_genes.gene_symbols),:)=0;
end
if size(scRNA_read_count_matrix,2)<length(scRNA_cells.cell_barcode)
    scRNA_read_count_matrix(:,end+1:length(scRNA_cells.cell_barcode))=0;
end
clear scRNA_read_count_table.peak_index scRNA_read_count_table.cell_index scRNA_read_count_table.readcount


%load scATACseq data
scATAC_assay_temp = readtable(scATAC_peak_file_path);
scATAC_peaks.peak_index=scATAC_assay_temp{:,1};
scATAC_peaks.chr=scATAC_assay_temp{:,2};
scATAC_peaks.point1=scATAC_assay_temp{:,3};
scATAC_peaks.point2=scATAC_assay_temp{:,4};
scATAC_peaks.chr_num=zeros(length(scATAC_peaks.chr), 1);
for i=1:22
    scATAC_peaks.chr_num(strcmp(scATAC_peaks.chr, ['chr', num2str(i)])>0,1)=i;
    %we exclude X and Y chromosome by setting its chr index as 0
end


[scATAC_cells.cell_index, scATAC_cells.cell_barcode, scATAC_cells.cell_type, scATAC_cells.subject_ID]=...
    textread(scATAC_cellmeta_file_path, '%d %s %s %s', 'delimiter', '\t');

[scATAC_read_count_table.peak_index, scATAC_read_count_table.cell_index, scATAC_read_count_table.readcount]=textread(scATAC_readcount_file_path, '%d %d %d');
scATAC_read_count_matrix=sparse(scATAC_read_count_table.peak_index,scATAC_read_count_table.cell_index,scATAC_read_count_table.readcount);
if size(scATAC_read_count_matrix,1)<length(scATAC_peaks.chr)
    scATAC_read_count_matrix(end+1:length(scATAC_peaks.chr),:)=0;
end
if size(scATAC_read_count_matrix,2)<length(scATAC_cells.cell_barcode)
    scATAC_read_count_matrix(:,end+1:length(scATAC_cells.cell_barcode))=0;
end
clear scATAC_read_count_table.peak_index scATAC_read_count_table.cell_index scATAC_read_count_table.readcount



%load TF motif prior
[Motifs.motif_index, Motifs.name]=textread(Motif_name_file_path, '%d %s');

[Motif_mapping_table.peak_index, Motif_mapping_table.motif_index, Motif_mapping_table.flag]=textread(Motif_mapping_file_path, '%d %d %d');
TF_peak_binding_matrix=sparse(Motif_mapping_table.peak_index, Motif_mapping_table.motif_index, Motif_mapping_table.flag);
if size(TF_peak_binding_matrix,1)<length(scATAC_peaks.chr)
    TF_peak_binding_matrix(end+1:length(scATAC_peaks.chr),:)=0;
end
if size(TF_peak_binding_matrix,2)<length(Motifs.name)
    TF_peak_binding_matrix(:,end+1:length(Motifs.name))=0;
end
clear Motif_mapping_table.peak_index Motif_mapping_table.motif_index Motif_mapping_table.flag


%load TAD prior
if TAD_flag==1%TAD_flag=1 if TAD prior is provided, TAD_flag=0 if not
    [TAD.chr, TAD.left_boundary, TAD.right_boundary]=textread(TAD_file_path, '%s %d %d');
    TAD.chr_num=zeros(length(TAD.chr), 1);
    for i=1:22
        TAD.chr_num(strcmp(TAD.chr, ['chr', num2str(i)])>0,1)=i;
        %we exclude X and Y chromosome by setting its chr index as 0
    end
else
    Distance_control=5e5;%500kb to TSS
end

[Refseq.chr, Refseq.strand, Refseq.start, Refseq.end, Refseq.gene_name]=textread(Ref_seq_file_path, '%s %s %d %d %s', 'headerlines', 1);
Refseq.chr_num=zeros(length(Refseq.chr), 1);
for i=1:22
    Refseq.chr_num(strcmp(Refseq.chr, ['chr', num2str(i)])>0,1)=i;
    %we exclude X and Y chromosome by setting its chr index as 0
end


%summary of loaded data
Group_conditions = unique(scRNA_cells.condition);
fprintf('We detected %d conditions from the meta file.\n\n', length(Group_conditions))

scRNA_samples = unique(scRNA_cells.subject_ID);
fprintf('The input scRNAseq data includes %d genes, %d cells from %d samples/subecjts.\n\n', length(scRNA_genes.gene_symbols), length(scRNA_cells.cell_barcode), length(scRNA_samples))

scATAC_samples = unique(scATAC_cells.subject_ID);
fprintf('The input scATACseq data includes %d peaks, %d cells from %d samples/subecjts.\n\n', length(scATAC_peaks.peak_index), length(scATAC_cells.cell_barcode), length(scATAC_samples))

Common_samples=intersect(scRNA_samples, scATAC_samples);
fprintf('There are paired data for %d samples/subecjts. (check sample IDs if this number is lower than expected)\n\n', length(Common_samples))


fprintf('%d motifs, %d candidate chromatin sites and %d candidate genes are provided.\n\n\n\n', length(Motifs.name), length(Candidate_peaks.chr), length(Candidate_genes.gene_symbols))



%********************** candidate circuits construction ****************


if TAD_flag==1
[Candidate_TFs, Candidate_TF_log2Count,...
    Candidate_peaks, Candidate_Peak_log2Count,...
    Candidate_genes, Candidate_Gene_log2Count,...
    Candidate_TF_Peak_Binding, Candidate_Peak_Gene_looping]=...
    Candidate_circuits_construction_with_TAD(Common_samples, Candidate_genes, Candidate_peaks,...
                                  scRNA_genes, scRNA_cells, scRNA_read_count_matrix, ...
                                  scATAC_peaks, scATAC_cells, scATAC_read_count_matrix,...
                                  Motifs, TF_peak_binding_matrix,...
                                  Refseq, TAD);
else
[Candidate_TFs, Candidate_TF_log2Count,...
    Candidate_peaks, Candidate_Peak_log2Count,...
    Candidate_genes, Candidate_Gene_log2Count,...
    Candidate_TF_Peak_Binding, Candidate_Peak_Gene_looping]=...
    Candidate_circuits_construction_without_TAD(Common_samples, Candidate_genes, Candidate_peaks,...
                                  scRNA_genes, scRNA_cells, scRNA_read_count_matrix, ...
                                  scATAC_peaks, scATAC_cells, scATAC_read_count_matrix,...
                                  Motifs, TF_peak_binding_matrix,...
                                  Refseq, Distance_control);
end


%********************** Model initialization ****************************


fprintf('MAGICAL model initialization ...\n\n')

S=length(Common_samples);
M=length(Candidate_TFs);
P=length(Candidate_peaks.peak_index);
G=length(Candidate_genes.gene_symbols);


[T_prior, T_mean, T_var, B_prior, B_mean, B_var, B_prob, L_prior, L_mean, L_var, L_prob]=...
    MAGICAL_initialization(Candidate_TF_log2Count, Candidate_Peak_log2Count, Candidate_Gene_log2Count, Candidate_TF_Peak_Binding, Candidate_Peak_Gene_looping, M);


%************************* MAGICAL sampling ****************************


fprintf('MAGICAL work starts ...\n\n');

A=Candidate_Peak_log2Count;
R=Candidate_Gene_log2Count;
B_state=full(Candidate_TF_Peak_Binding);
L_state=full(Candidate_Peak_Gene_looping);

[Candidate_TF_Peak_Binding_prob, Candidate_Peak_Gene_Looping_prob]=...
    MAGICAL_estimation(A, R, B_state, L_state, T_prior, T_mean, T_var, B_prior, B_mean, B_var, B_prob, L_prior, L_mean, L_var, L_prob, S, P, G, iteration_num);



%************************* MAGICAL output ******************************


fid=fopen(Output_file_path, 'w');
fprintf(fid, 'Gene_symbol\tGene_chr\tGene_TSS\tPeak_chr\tPeak_start\tPeak_end\tLooping_prob\tTFs(binding prob)\n');

[xx,yy]=find(Candidate_Peak_Gene_Looping_prob>prob_threshold_peak_gene_looping);

TF_vector=zeros(length(Candidate_TFs), 1);
for i=1:length(xx)
    [TF_prob, TF_index]=sort(full(Candidate_TF_Peak_Binding_prob(xx(i),:)), 'descend');
    index=find(TF_prob>prob_threshold_TF_peak_binding);
    if ~isempty(index)
        
        fprintf(fid, '%s\tchr%d\t%d\tchr%d\t%d\t%d\t%G\t', ...
            Candidate_genes.gene_symbols{yy(i)}, Candidate_genes.gene_TSS(yy(i),:),...
            Candidate_peaks.chr_num(xx(i)), Candidate_peaks.point1(xx(i)),  Candidate_peaks.point2(xx(i)), full(Candidate_Peak_Gene_Looping_prob(xx(i),yy(i))));
        
        for t=1:length(index)
            fprintf(fid, '%s (%.2f), ', Candidate_TFs{TF_index(t)}, TF_prob(t));
            TF_vector(TF_index(t))=TF_vector(TF_index(t))+1;
        end
        
        fprintf(fid, '\n');
    end
end
fprintf('MAGICAL selected regulatory circuits with %d TFs, %d peaks and %d genes.\n\n\n\n', length(find(TF_vector>1)), length(unique(xx)), length(unique(yy)))





