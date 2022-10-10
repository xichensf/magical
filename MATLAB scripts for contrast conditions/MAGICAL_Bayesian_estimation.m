function MAGICAL_post = MAGICAL_Bayesian_estimation(Candidate_TFs, Candidate_Peaks, Candidate_Genes, Candidate_Gene_TSS, ...
    Candidate_Peak_log2Count, Candidate_Gene_log2Count, Candidate_TF_Peak_Binding, Candidate_Peak_Gene_looping, ...
                            P_prior, P_mean, P_var, B_prior, B_mean, B_var, B_prob, L_prior, L_mean, L_var, L_prob, S, T, F, G)

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
    iteration_num=1000;
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
        
        %Step 8: Sample Summary
        B_state_frq=B_state_frq+B_state;
        L_state_frq=L_state_frq+L_state;
        
        if mod(i,iteration_seg)==0
            fprintf(2, 'MAGICAL finished %d percent\n\n', 10*i/iteration_seg)
        end
    end
    
    
    % output posterior TFs, Peaks, Genes, TF-Peak binding and Peak-Gene looping
    
    MAGICAL_post.TFs=Candidate_TFs;
    MAGICAL_post.Peaks=Candidate_Peaks;
    MAGICAL_post.Genes=Candidate_Genes;
    MAGICAL_post.Gene_TSS=Candidate_Gene_TSS;
    MAGICAL_post.TF_Peak_Binding_prob=B_prob;
    MAGICAL_post.Peak_Gene_Looping_prob=L_prob;
