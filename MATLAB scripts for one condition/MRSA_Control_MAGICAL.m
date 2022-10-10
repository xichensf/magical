% /Applications/MATLAB_R2019a.app/bin/matlab -batch

clear
close all hidden
clc
warning('off','all')

% very important. please chnage chr1 - chr22 to 1-22 for quick indexing,also chnage chrX to 23 and chr to 24
[Peaks(:,1), Peaks(:,2), Peaks(:,3)]=textread('Peaks.txt', '%d %d %d', 'headerlines', 1);

ATAC_Cell_types={'B memory', 'B naive', 'CD4 Naive', 'CD4 TCM', 'CD8 Naive', 'CD8 TEM', 'CD14 Mono', 'CD16 Mono', 'cDC2', 'MAIT', 'NK', 'NK_CD56bright', 'Platelet'};
ATAC_Cell_type_proportion = readtable('ATAC_cell_type_proportion_20000peak_TSS12_Nu2_reso2_L2_manual.csv');
Sample_meta.ATAC_sample_names = ATAC_Cell_type_proportion.Sample_name;
Sample_meta.ATAC_sample_condition = ATAC_Cell_type_proportion.Condition;


% the order of RNA cell types must be the same as ATAC cell types
RNA_Cell_types ={'B_memory', 'B_naive', 'CD4_Naive', 'CD4_TCM', 'CD8_Naive', 'CD8_TEM_CTL', 'CD14_Mono', 'CD16_Mono', 'cDC2', 'MAIT', 'NK', 'NK_CD56bright', 'Platelet'};
RNA_Cell_type_proportion = readtable('RNA_MSSA-MRSA-CTRL_cell-type-proportions_wide_070621.csv');
Sample_meta.RNA_sample_names = RNA_Cell_type_proportion.Sample_name;
Sample_meta.RNA_sample_condition = RNA_Cell_type_proportion.Condition;


% load TF motif-peak mapping
peakmotifmatches=readtable('peak_motif_matches.txt');
TFs=peakmotifmatches.Properties.VariableNames(4:end);
TF_peaks=peakmotifmatches{:,1:3};
TF_Peak_Binding=peakmotifmatches{:,4:end};
clear peakmotifmatches


%MRSA vs Control 
fid=fopen('MRSA_Control_MAGICAL_temp.txt', 'w');
fprintf(fid, 'Cell_type\tGene_symbol\tGene_chr\tGene_TSS\tGene_Case_pct\tGGene_Control_pct\tene_sc_Log2FC\tGene_sc_p_val_adj\tGene_bulk_Log2FC\tGene_bulk_pval\tPeak_chr\tPeak_start\tPeak_end\tPeak_sc_log2FC\tPeak_sc_p_val_adj\tPeak_bulk_Log2FC\tPeak_bulk_pval\tATAC_RNA_Bayesian_prob\tTFs\n');


for c=1:length(ATAC_Cell_types)
    fprintf('Integrate ATAC and RNA data for %s:\n\n', ATAC_Cell_types{c})
    
    %load peak scATAC differential statistics
    Peak_sc_FDR = readmatrix(['scATAC_DAPs/', ATAC_Cell_types{c},'/MRSA_Control_diff.xlsx'],'Sheet','FDR');
    Peak_sc_Log2FC= readmatrix(['scATAC_DAPs/', ATAC_Cell_types{c},'/MRSA_Control_diff.xlsx'],'Sheet','Log2FC');
    if length(Peak_sc_FDR)<length(Peak_sc_Log2FC)
        Peak_sc_FDR(end+1:length(Peak_sc_Log2FC))=1;
    end     
    Diff_peak_index=find(abs(Peak_sc_Log2FC)>0.1 & Peak_sc_FDR<0.05);
    
    sc_ATAC.Peaks=Peaks(Diff_peak_index,:); 
    sc_ATAC.FDR=Peak_sc_FDR(Diff_peak_index); 
    sc_ATAC.Log2FC=Peak_sc_Log2FC(Diff_peak_index);
    
    Diff_Peaks=sc_ATAC.Peaks;
    
    fprintf(2, '%d DA peaks from single cell differential analysis\n\n', length(Diff_Peaks))
    
    
    
     %% cell type-specific pseudo bulk ATAC
    fprintf('ATAC psuedo-bulk differential analysis\n\n')
    ATAC_bulk = readtable(['scATAC_pseudobulk/pseudo_bulk_ATAC_count_', ATAC_Cell_types{c}, '.txt']);
    ATAC_samples = ATAC_bulk.Properties.VariableNames(2:end);
    [ATAC_Case_samples,ATAC_Case_sample_idx]=intersect(ATAC_samples, Sample_meta.ATAC_sample_names(strcmp(Sample_meta.ATAC_sample_condition, 'MRSA')>0), 'stable');
    [ATAC_Control_samples,ATAC_Control_sample_idx]=intersect(ATAC_samples, Sample_meta.ATAC_sample_names(strcmp(Sample_meta.ATAC_sample_condition, 'Control')>0), 'stable');
    
    [Pseudo_ATAC.Peaks, Pseudo_ATAC.log2FC, Pseudo_ATAC.pvalue, Pseudo_ATAC.robust_pvalue] = ATAC_psuedo_bulk_differential(Diff_Peaks, Peaks, ATAC_bulk, ATAC_Case_sample_idx', ATAC_Control_sample_idx');
    Diff_Peaks=Pseudo_ATAC.Peaks(abs(Pseudo_ATAC.log2FC)>0.25 & (Pseudo_ATAC.pvalue<0.05 | Pseudo_ATAC.robust_pvalue<0.05),:);
    fprintf(2, '%d DA peaks are still signficant in psuedo bulk\n\n', length(Diff_Peaks))    
    
    
    
    
    %% Enrichment analysis on differential peaks (with active peaks as background), to select candidate TFs regulating these peaks
    fprintf('Select candidate TFs that have motifs mapped in differential peaks\n\n')
    fprintf(2, 'scanned %d TF motifs at peak regions\n\n', length(TFs))
  
    
    
    
    %% scRNA differential statistics
    fprintf('load scRNA differential statistics from R output\n\n')
    [DEG.symbol, DEG.p_val, DEG.log2FC, DEG.pct_1, DEG.pct_2, DEG.p_val_adj]=...
        textread(['scRNA_DEGs/MRSA_vs_CTRL/MRSA-vs-CTRL-degs-', RNA_Cell_types{c}, '-0627.csv'], '%s %f %f %f %f %f', 'delimiter', ',', 'headerlines', 1);
    
    load Protein_coding_genes.mat
    [aa, bidex]=intersect(DEG.symbol, Gene_protein_coding, 'stable');
    DEG.symbol=DEG.symbol(bidex);
    DEG.p_val=DEG.p_val(bidex);
    DEG.log2FC=DEG.log2FC(bidex); 
    DEG.pct_1=DEG.pct_1(bidex); 
    DEG.pct_2=DEG.pct_2(bidex); 
    DEG.p_val_adj=DEG.p_val_adj(bidex);
    
    
    Diff_gene_index=find(DEG.p_val_adj<0.05 & abs(DEG.log2FC)>0.1 & max([DEG.pct_1, DEG.pct_2], [], 2)>0.1);
    
    sc_RNA.Genes=DEG.symbol(Diff_gene_index);
    sc_RNA.pct_1=DEG.pct_1(Diff_gene_index); 
    sc_RNA.pct_2=DEG.pct_2(Diff_gene_index);
    sc_RNA.avg_log2FC=DEG.log2FC(Diff_gene_index);
    sc_RNA.p_val=DEG.p_val(Diff_gene_index);
    sc_RNA.p_val_adj=DEG.p_val_adj(Diff_gene_index);
 
    Diff_Genes = sc_RNA.Genes;
    
    fprintf(2, '%d DE genes from single cell differential analysis\n\n', length(Diff_Genes))
    clear DEG
    
    
    
    
    %% pseudo bulk RNA differential statistics        
    fprintf('RNA psuedo-bulk differential analysis\n\n')
    RNA_bulk = readtable(['scRNA_pseudobulk/per-celltype/MRSA-MSSA-pseudobulk-per-samples-', RNA_Cell_types{c}, '-0624.csv']);

    RNA_samples = RNA_bulk.Properties.VariableNames(2:end);
    
    [RNA_Case_samples,RNA_Case_sample_idx]=intersect(RNA_samples, Sample_meta.RNA_sample_names(strcmp(Sample_meta.RNA_sample_condition, 'MRSA')>0), 'stable');
    [RNA_Control_samples,RNA_Control_sample_idx]=intersect(RNA_samples, Sample_meta.RNA_sample_names(strcmp(Sample_meta.RNA_sample_condition, 'Control')>0), 'stable');
    
    [Pseudo_RNA.Genes, Pseudo_RNA.log2FC, Pseudo_RNA.pvalue, Pseudo_RNA.robust_pvalue] = RNA_psuedo_bulk_differential(Diff_Genes, RNA_bulk, RNA_Case_sample_idx', RNA_Control_sample_idx');
    Diff_Genes=Pseudo_RNA.Genes(abs(Pseudo_RNA.log2FC)>0.25 & (Pseudo_RNA.pvalue<0.05 | Pseudo_RNA.robust_pvalue<0.05));
    fprintf(2, '%d DE genes are still significant in psuedo bulk\n\n', length(Diff_Genes))
    

    
    
    %% MAGICAL integration
    fprintf('MAGICAL integration starts!\n\n')
    %     MAGICAL_post=MAGICAL_V1(Candidate_Peak_ID, Candidate_Peaks, ATAC, Candidate_Gene_Symbol, RNA, Candidate_TFs, Candidate_TF_Peak_Binding, Sample_meta);
    %************************* Input data *********************
    % Candidate_Peak_ID: P x 1 string vector
    % Candidate_Peaks: P x 3 matrix with chr, point1, point2
    % ATAC: P x S matrix, raw pseudo bulk ATAC count for all peaks and S samples
    % Candidate_Gene_Symbol: G x 1 string vector with gene names
    % RNA: G x S matrix, raw pseudo bulk RNA count for all genes and S samples
    % Candidate_TFs: T x 1 string vector for enriched TFs
    % Candidate_TF_Peak_Binding: P x T binary matrix, binding state for P peaks and T TFs
    
    %********** Initial Integration to select TFs, peaks and genes ***********
    [ATAC_Case_samples,ATAC_Case_sample_idx]=intersect(ATAC_samples, Sample_meta.ATAC_sample_names(strcmp(Sample_meta.ATAC_sample_condition, 'Control')==0), 'stable');
    [ATAC_Control_samples,ATAC_Control_sample_idx]=intersect(ATAC_samples, Sample_meta.ATAC_sample_names(strcmp(Sample_meta.ATAC_sample_condition, 'Control')>0), 'stable');
    
    [RNA_Case_samples,RNA_Case_sample_idx]=intersect(RNA_samples, Sample_meta.RNA_sample_names(strcmp(Sample_meta.RNA_sample_condition, 'Control')==0), 'stable');
    [RNA_Control_samples,RNA_Control_sample_idx]=intersect(RNA_samples, Sample_meta.RNA_sample_names(strcmp(Sample_meta.RNA_sample_condition, 'Control')>0), 'stable');

    [Common_samples, Candidate_TFs, Candidate_TF_log2Count,Candidate_Peaks, Candidate_Peak_log2Count,...
        Candidate_Genes, Candidate_Gene_TSS, Candidate_Gene_log2Count,Candidate_TF_Peak_Binding, Candidate_Peak_Gene_looping] =...
                   Initial_integration(Diff_Peaks, Peaks, ATAC_bulk, ATAC_Case_samples, ATAC_Control_samples, Diff_Genes, RNA_bulk, RNA_Case_samples, RNA_Control_samples,TFs, TF_peaks, TF_Peak_Binding);
    
    S=length(Common_samples);
    T=length(Candidate_TFs);
    F=length(Candidate_Peaks);
    G=length(Candidate_Genes);
    
    fprintf(2, 'Initial integration associated %d TFs, %d DA peaks, and %d DE genes together\n\n', T, F, G)
    
    if G==0
        fprintf(2, 'Too few candidate peaks and genes (<5), MAGICAL skipped the current cell type!\n\n\n\n')
        continue;
    else
        
        
    %********** Model variable prior calculation *****************************
    %Note: in MCMC, initial values may not be that important, but using more
    %importantive prior could speed up the sampling convergence
    
    fprintf('MAGICAL model initialization\n\n')
    
    [P_prior, P_mean, P_var, B_prior, B_mean, B_var, B_prob, L_prior, L_mean, L_var, L_prob] = ...
        MAGICAL_prior(Candidate_TF_log2Count,Candidate_Peak_log2Count,Candidate_Gene_log2Count,Candidate_TF_Peak_Binding, Candidate_Peak_Gene_looping, S, T, F, G);

    
   %********** MAGICAL integration running here *****************************

    MAGICAL_post = MAGICAL_Bayesian_estimation(Candidate_TFs, Candidate_Peaks, Candidate_Genes, Candidate_Gene_TSS, ...
                            Candidate_Peak_log2Count, Candidate_Gene_log2Count, Candidate_TF_Peak_Binding, Candidate_Peak_Gene_looping, ...
                            P_prior, P_mean, P_var, B_prior, B_mean, B_var, B_prob, L_prior, L_mean, L_var, L_prob, S, T, F, G);


   %********** write results to output files *****************************
                        
    [xx,yy]=find(MAGICAL_post.Peak_Gene_Looping_prob>0.95);
    
    for i=1:length(xx)
        Gene_index_1=find(strcmp(sc_RNA.Genes, MAGICAL_post.Genes(yy(i)))>0);
        Gene_index_2=find(strcmp(Pseudo_RNA.Genes, MAGICAL_post.Genes(yy(i)))>0);
        
        [aa, Peak_index_1]=intersect(sc_ATAC.Peaks, MAGICAL_post.Peaks(xx(i), :), 'rows');
        [bb, Peak_index_2]=intersect(Pseudo_ATAC.Peaks, MAGICAL_post.Peaks(xx(i),:), 'rows');
        
        fprintf(fid, '%s\t%s\t%d\t%d\t%f\t%f\t%f\t%G\t%f\t%G\t%d\t%d\t%d\t%f\t%G\t%f\t%G\t%G\t', ...
            ATAC_Cell_types{c}, MAGICAL_post.Genes{yy(i)}, MAGICAL_post.Gene_TSS(yy(i), :),...
            sc_RNA.pct_1(Gene_index_1(1)), sc_RNA.pct_2(Gene_index_1(1)),...
            sc_RNA.avg_log2FC(Gene_index_1(1)), sc_RNA.p_val_adj(Gene_index_1(1)),...
            Pseudo_RNA.log2FC(Gene_index_2(1)), min(Pseudo_RNA.pvalue(Gene_index_2(1)), Pseudo_RNA.robust_pvalue(Gene_index_2(1))),...
            MAGICAL_post.Peaks(xx(i),:),...
            sc_ATAC.Log2FC(Peak_index_1(1)), sc_ATAC.FDR(Peak_index_1(1)),...
            Pseudo_ATAC.log2FC(Peak_index_2(1)), min(Pseudo_ATAC.pvalue(Peak_index_2(1)), Pseudo_ATAC.robust_pvalue(Peak_index_2(1))),...
            full(MAGICAL_post.Peak_Gene_Looping_prob(xx(i),yy(i))));
        
        [TF_prob, TF_index]=sort(full(MAGICAL_post.TF_Peak_Binding_prob(xx(i),:)), 'descend');
        for t=1:length(TF_index)
            if TF_prob(t)>0.8
                fprintf(fid, '%s(%.2f), ', MAGICAL_post.TFs{TF_index(t)}, TF_prob(t));
            end
        end
        fprintf(fid, '\n');
    end
    fprintf('MAGICAL finally integrated %d peaks and %d genes for %s .\n\n\n\n', length(unique(xx)), length(unique(yy)), ATAC_Cell_types{c})
    end
end
    
