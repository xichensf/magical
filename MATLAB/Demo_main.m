% /Applications/MATLAB_R2019a.app/bin/matlab -batch Demo_main

clear
close all hidden
clc
warning('off','all')
%load candidate genes






%load candidate peaks


[Candidate_Peaks(:,1), Candidate_Peaks(:,2), Candidate_Peaks(:,3)]=textread('Peaks.txt', '%s %d %d', 'headerlines', 1);



%load scRNAseq data






%load scATACseq data


[scATAC_Peaks(:,1), scATAC_Peaks(:,2), scATAC_Peaks(:,3)]=textread('Peaks.txt', '%s %d %d', 'headerlines', 1);
Sample_meta = readtable(' Metadata for BSU Samples 6-3-21_TC.txt');



%load TF-peak binding prior







%load TAD prior or distance to TSS control: TAD_flag=1 or 0

TAD_flag=1;










for c=3:3%size(Cell_type_file_list)
    fprintf('Integrate ATAC and RNA data for %s:\n\n', Cell_type_file_list{c,1})
    
    % select candidate ATAC peaks for case vs control
    fprintf('load scATAC differential statistics from R output\n\n')
    %connect to R output
    Peak_sc_diff = readtable(['scATAC/Differential_peaks/Exercise_post_pre/', Cell_type_file_list{c,2}]);
    index=find(strcmp(Peak_sc_diff.cluster, 'Pre')>0);
    Peak_sc_diff.avg_log2FC(index)=-1*Peak_sc_diff.avg_log2FC(index);
    DA_peak_sc_index=find(abs(Peak_sc_diff.avg_log2FC)>0.1 & Peak_sc_diff.p_val_adj<0.05 & (Peak_sc_diff.pct_1>0.1 | Peak_sc_diff.pct_2>0.1));
    Peak_diff_ID=Peak_sc_diff.Var1(DA_peak_sc_index);
    Active_peak_sc_index=find((Peak_sc_diff.pct_1>0.1 | Peak_sc_diff.pct_2>0.1));
    Peak_active_ID=Peak_sc_diff.Var1(Active_peak_sc_index);
    
    Peak_control_diff = readtable(['scATAC/Differential_peaks/Control_post_pre/', Control_file_list{c,1}]);
    Background_peak_index=find(abs(Peak_control_diff.avg_log2FC)>0.1 & Peak_control_diff.p_val_adj<0.05 &...
        (Peak_control_diff.pct_1>0.1 | Peak_control_diff.pct_2>0.1));
    Background_sc_diff_ID=Peak_control_diff.Var1(Background_peak_index);
    Peak_diff_ID=setdiff(Peak_diff_ID, Background_sc_diff_ID);
    
    fprintf(2, '%d DA peaks from single cell differential analysis\n\n', length(Peak_diff_ID))
    
    
    %% pseudo bulk ATAC differential statistics
    fprintf('ATAC psuedo-bulk differential analysis\n\n')
    ATAC = readtable(['scATAC/Pseudo_bulk/', Cell_type_file_list{c,3}]);
    ATAC_samples = ATAC.Properties.VariableNames(2:end);
    [ATAC_Case_samples,ATAC_Case_sample_idx]=intersect(ATAC_samples, Sample_meta.biospecimenID(strcmp(Sample_meta.Group, 'Exercise')>0), 'stable');
    [ATAC_Control_samples,ATAC_Control_sample_idx]=intersect(ATAC_samples, Sample_meta.biospecimenID(strcmp(Sample_meta.Group, 'Baseline')>0), 'stable');
    
    [ATAC_Peak_ID, ATAC_log2FC, ATAC_pvalue, ATAC_robust_pvalue] = ATAC_psuedo_bulk_differential(Peak_diff_ID, ATAC, ATAC_Case_sample_idx, ATAC_Control_sample_idx);
    DA_peak_bulk_index=find(abs(ATAC_log2FC)>0.3 & (ATAC_pvalue<0.05 | ATAC_robust_pvalue<0.05));
    Peak_diff_ID=ATAC_Peak_ID(DA_peak_bulk_index);
    fprintf(2, '%d DA peaks are still signfiicant in psuedo bulk\n\n', length(Peak_diff_ID))
    
    
    %% Enrichment analysis on differential peaks (with active peaks as background), to select candidate TFs regulating these peaks
    fprintf('Select candidate TFs that are enriched in differential peaks, with active peaks as background\n\n')
    TF_binding=load('Muscle_Peak_Motifs_mapping');
    [Candidate_Peak_ID,bidex]=intersect(TF_binding.Peak_ID, Peak_diff_ID);
    Motif_diff_peak_binding=TF_binding.Peak_Motif_binding(bidex,:);
    
    [Peak_active_ID,bidex]=intersect(TF_binding.Peak_ID, Peak_active_ID);
    Motif_active_peak_binding=TF_binding.Peak_Motif_binding(bidex,:);
    
    
    TF_num=full(sum(Motif_diff_peak_binding));
    TF_pct=full(sum(Motif_diff_peak_binding))/length(Candidate_Peak_ID);
    TF_enrichment_pvalue=ones(1,length(length(TF_binding.Motifs)));
    TF_enrichment_FC=zeros(1,length(length(TF_binding.Motifs)));
    for t=1:length(TF_binding.Motifs)
        TF_enrichment_FC(t)=(TF_num(t)/length(Candidate_Peak_ID))/(sum(Motif_active_peak_binding(:,t))/length(Peak_active_ID));
    end
    TF_index=find(TF_pct>0.1 & TF_num>30 & TF_enrichment_FC>1.5);
    
    if ~isempty(TF_index)
        Candidate_TFs=TF_binding.Motifs(TF_index);
        Candidate_TF_Peak_Binding=Motif_diff_peak_binding(:,TF_index);
    else
        fprintf('Too few peaks with TF binding sites.\nMAGICAL not applicable to this cell type!\n\n\n\n')
%         continue;
    end
    
    Binding_peak_index=find(sum(Candidate_TF_Peak_Binding, 2)>0);
    Candidate_Peak_ID=Candidate_Peak_ID(Binding_peak_index);
    Candidate_TF_Peak_Binding=Candidate_TF_Peak_Binding(Binding_peak_index,:);
    
    [Candidate_Peak_ID,bidex,cidex]=intersect(Candidate_Peak_ID, Peak_ID, 'stable');
    Candidate_TF_Peak_Binding=Candidate_TF_Peak_Binding(bidex,:);
    Candidate_Peaks=Peaks(cidex,:);
    fprintf(2, '%d enriched TFs\n\n', length(Candidate_TFs))
    
    
    
    %% scRNA differential statistics
    fprintf('load scRNA differential statistics from R output\n\n')
    RNA_sc_diff = readtable(['scRNA/Differential_genes/Exercise_post_pre/', Cell_type_file_list{c,4}]);
    index=find(strcmp(RNA_sc_diff.cluster, 'Pre')>0);
    RNA_sc_diff.avg_log2FC(index)=-1*RNA_sc_diff.avg_log2FC(index);
    DE_gene_sc_index=find(RNA_sc_diff.p_val_adj<0.05 & abs(RNA_sc_diff.avg_log2FC)>0.1 & max([RNA_sc_diff.pct_1, RNA_sc_diff.pct_2], [], 2)>0.1);
    Gene_diff_symbol=RNA_sc_diff.Var1(DE_gene_sc_index);
    
    Gene_control_diff = readtable(['scRNA/Differential_genes/Control_post_pre/', Control_file_list{c,2}]);
    Background_gene_index=find(Gene_control_diff.p_val_adj<0.05 & abs(Gene_control_diff.avg_log2FC)>0.1 &...
        max([Gene_control_diff.pct_1, Gene_control_diff.pct_2], [], 2)>0.1);
    Background_sc_diff_symbol=Gene_control_diff.Var1(Background_gene_index);
    Gene_diff_symbol=setdiff(Gene_diff_symbol, Background_sc_diff_symbol);
    
    fprintf(2, '%d DE genes from single cell differential analysis\n\n', length(Gene_diff_symbol))
    
    
    %% pseudo bulk RNA differential statistics
    fprintf('RNA psuedo-bulk differential analysis\n\n')
    RNA = readtable(['scRNA/Pseudo_bulk/', Cell_type_file_list{c,5}]);
    RNA_samples = RNA.Properties.VariableNames(2:end);
    [RNA_Case_samples,RNA_Case_sample_idx]=intersect(RNA_samples, Sample_meta.biospecimenID(strcmp(Sample_meta.Group, 'Exercise')>0), 'stable');
    [RNA_Control_samples,RNA_Control_sample_idx]=intersect(RNA_samples, Sample_meta.biospecimenID(strcmp(Sample_meta.Group, 'Baseline')>0), 'stable');
    
    [Gene_names, Gene_log2FC, Gene_pvalue, Gene_robust_pvalue]=RNA_psuedo_bulk_differential(Gene_diff_symbol, RNA, RNA_Case_sample_idx, RNA_Control_sample_idx);
    DE_gene_bulk_index=find(abs(Gene_log2FC)>0.5 & (Gene_pvalue<0.05 | Gene_robust_pvalue<0.05));
    Candidate_Gene_Symbol=Gene_names(DE_gene_bulk_index);
    fprintf(2, '%d DE genes are still significant in psuedo bulk\n\n', length(Candidate_Gene_Symbol))
    
    
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
    [Common_samples, Candidate_TFs, Candidate_TF_log2Count,...
        Candidate_Peak_ID, Candidate_Peaks, Candidate_Peak_log2Count,...
        Candidate_Gene_Symbol, Candidate_Gene_TSS, Candidate_Gene_log2Count,...
        Candidate_TF_Peak_Binding, Candidate_Peak_Gene_looping]= Initial_integration(Candidate_Peak_ID, Candidate_Peaks, ATAC,...
        Candidate_Gene_Symbol, RNA, Candidate_TFs, Candidate_TF_Peak_Binding, Sample_meta);
    
    S=length(Common_samples);
    T=length(Candidate_TFs);
    F=length(Candidate_Peaks);
    G=length(Candidate_Gene_Symbol);
    
    fprintf(2, 'Initial integration associated %d TFs, %d DA peaks, and %d DE genes together\n\n', T, F, G)
    
    %********** Model variable prior calculation *****************************
    %Note: in MCMC, initial values may not be that important, but using more
    %importantive prior could speed up the sampling convergence
    
    fprintf('MAGICAL model initialization\n\n')
    
    [P_prior, P_mean, P_var, B_prior, B_mean, B_var, B_prob, L_prior, L_mean, L_var, L_prob]=...
        MAGICAL_prior(Candidate_TF_log2Count,Candidate_Peak_log2Count,Candidate_Gene_log2Count,Candidate_TF_Peak_Binding, Candidate_Peak_Gene_looping, S, T, F, G);
    
    
    %********** Initial round of sampling, using their prior values **********
    
    A=Candidate_Peak_log2Count;
    R=Candidate_Gene_log2Count;
    
    P=P_prior;
    B=B_prior;
    L=L_prior;
    
    B_state=full(Candidate_TF_Peak_Binding);
    L_state=full(Candidate_Peak_Gene_looping);
    
    B_state_frq=full(Candidate_TF_Peak_Binding);
    L_state_frq=full(Candidate_Peak_Gene_looping);
    
    alpha_A=1;
    beta_A=1;
    ATAC_fitting_residue=A-B*P;
    sigma_A_noise=var(ATAC_fitting_residue(:));
    
    % RNA_fitting_residue=R-L'*(B*P);
    % sigma_R_noise=var(RNA_fitting_residue(:));
    alpha_R=1;
    beta_R=1;
    RNA_fitting_residue=R-L'*A;
    sigma_R_noise=var(RNA_fitting_residue(:));
    
    %************************* MAGICAL sampling ******************************
    iteration_num=10000;
    iteration_seg=iteration_num/10;
    for i=1:iteration_num
        
        %********** scATAC-seq fitting and variable sampling  ****************
        
        %Step 1: TF activity sampling
        P = TF_activity_P_sampling(A, B, P, P_mean, P_var, sigma_A_noise); 
        
        %Step 2: TF-peak binding weight sampling
        B = TF_peak_binding_B_sampling(A, B, P, B_state, B_mean, B_var, sigma_A_noise);
        
        %Step 3: TF-peak binding state update
        [B_state, B] = TF_peak_binary_binding_B_state_sampling(A, B, P, B_state, B_mean, B_var, B_prob, sigma_A_noise);
        
        %Step 4: ATAC fitting residue variance control
        ATAC_fitting_residue=A-B*P;
        sigma_A_noise = 1/gamrnd(alpha_A+1/2,1/(beta_A+sum(sum(ATAC_fitting_residue.^2))/(2*F*S)));
%         sigma_A_noise = (beta_A+sum(sum(ATAC_fitting_residue.^2))/(2*F*S))/chi2rnd(2*alpha_A+1);
        aa(i)=sigma_A_noise;
        
        %********** scRNA-seq fitting and variable sampling  *****************
        A_estimate=B*P;% or true A
        %Step 5: Peak-Gene looping weight sampling
        L = Peak_gene_looping_L_samping(R, L, A_estimate, L_state, L_mean, L_var, sigma_R_noise);
        
        %Step 6: Peak-Gene looping state update
        [L_state, L]=Peak_gene_binary_looping_L_state_samping(R, L, A_estimate, L_state, L_mean, L_var, L_prob, sigma_R_noise);
        
        %Step 7: RNA fitting residue variance control
        RNA_fitting_residue=R-L'*A_estimate;
        sigma_R_noise = 1/gamrnd(alpha_R+1/2,1/(beta_R+sum(sum(RNA_fitting_residue.^2))/(2*G*S)));
%         sigma_R_noise = (beta_R+sum(sum(RNA_fitting_residue.^2))/(2*G*S))/chi2rnd(2*alpha_R+1);
        bb(i) = sigma_R_noise;
        
        %Step 8: Sample Summary
        B_state_frq=B_state_frq+B_state;
        L_state_frq=L_state_frq+L_state;
        
        if mod(i,iteration_seg)==0
            fprintf(2, 'MAGICAL finished %d percent\n\n', 10*i/iteration_seg)
        end
    end
    
    
    % output posterior TFs, Peaks, Genes, TF-Peak binding and Peak-Gene looping
    
    MAGICAL_post.TFs=Candidate_TFs;
    MAGICAL_post.Peak_ID=Candidate_Peak_ID;
    MAGICAL_post.Peaks=Candidate_Peaks;
    MAGICAL_post.Genes=Candidate_Gene_Symbol;
    MAGICAL_post.Gene_TSS=Candidate_Gene_TSS;
    MAGICAL_post.TF_Peak_Binding=B_prob;
    MAGICAL_post.Peak_Gene_Looping=L_prob;
    MAGICAL_post.TF_Peak_Binding_prob=B_state_frq/iteration_num;
    MAGICAL_post.Peak_Gene_Looping_prob=L_state_frq/iteration_num;
        
    fid=fopen([Cell_type_file_list{c,1},'_MAGICAL_',num2str(iteration_num) ,'_sampling.txt'], 'w');
    fprintf(fid, 'Cell_type\tGene_symbol\tGene_chr\tGene_TSS\tGene_pct_post\tGene_pct_pre\tGene_sc_Log2FC\tGene_sc_p_val_adj\tGene_bulk_Log2FC\tGene_bulk_pval\tPeak_ID\tPeak_chr\tPeak_start\tPeak_end\tPeak_sc_post\tPeak_sc_pre\tPeak_sc_log2FC\tPeak_sc_p_val_adj\tPeak_bulk_Log2FC\tPeak_bulk_pval\tATAC_RNA_Bayesian_pval\tTFs\n');
    
    [xx,yy]=find(MAGICAL_post.Peak_Gene_Looping_prob>0.8);
    
    for i=1:length(xx)
        Gene_index_1=find(strcmp(RNA_sc_diff.Var1, MAGICAL_post.Genes(yy(i)))>0);
        Gene_index_2=find(strcmp(Gene_names, MAGICAL_post.Genes(yy(i)))>0);
        
        Peak_index_1=find(strcmp(Peak_sc_diff.Var1, MAGICAL_post.Peak_ID(xx(i)))>0);
        Peak_index_2=find(strcmp(ATAC_Peak_ID, MAGICAL_post.Peak_ID(xx(i)))>0);
        
        fprintf(fid, '%s\t%s\tchr%d\t%d\t%f\t%f\t%f\t%G\t%f\t%G\t%s\tchr%d\t%d\t%d\t%f\t%f\t%f\t%G\t%f\t%G\t%G\t', ...
            Cell_type_file_list{c,1}, MAGICAL_post.Genes{yy(i)}, MAGICAL_post.Gene_TSS(yy(i), :),...
            RNA_sc_diff.pct_1(Gene_index_1(1)), RNA_sc_diff.pct_2(Gene_index_1(1)),...
            RNA_sc_diff.avg_log2FC(Gene_index_1(1)), RNA_sc_diff.p_val_adj(Gene_index_1(1)),...
            Gene_log2FC(Gene_index_2(1)), min(Gene_pvalue(Gene_index_2(1)), Gene_robust_pvalue(Gene_index_2(1))),...
            MAGICAL_post.Peak_ID{xx(i)}, MAGICAL_post.Peaks(xx(i),:),...
            Peak_sc_diff.pct_1(Peak_index_1(1)), Peak_sc_diff.pct_2(Peak_index_1(1)),...
            Peak_sc_diff.avg_log2FC(Peak_index_1(1)), Peak_sc_diff.p_val_adj(Peak_index_1(1)),...
            ATAC_log2FC(Peak_index_2(1)), min(ATAC_pvalue(Peak_index_2(1)), ATAC_robust_pvalue(Peak_index_2(1))),...
            full(MAGICAL_post.Peak_Gene_Looping_prob(xx(i),yy(i))));
        
        [TF_prob, TF_index]=sort(full(MAGICAL_post.TF_Peak_Binding_prob(xx(i),:)), 'descend');
        for t=1:length(TF_index)
            if TF_prob(t)>0.8
                fprintf(fid, '%s(%.2f), ', MAGICAL_post.TFs{TF_index(t)}, TF_prob(t));
            end
        end
        fprintf(fid, '\n');
    end
    fprintf('MAGICAL integrated %d peaks and %d genes together.\n\n', length(unique(xx)), length(unique(yy)))
    fprintf('Results in file: %s_MAGICAL_results.txt\n\n\n\n', Cell_type_file_list{c,1})
    
end
    
