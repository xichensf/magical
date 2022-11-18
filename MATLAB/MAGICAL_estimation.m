function [Candidate_TF_Peak_Binding_prob, Candidate_Peak_Gene_Looping_prob]=...
    MAGICAL_estimation(scATAC_read_count_matrix, ATAC_cell_vector, Candidate_Peak_log2Count, scRNA_read_count_matrix, RNA_cell_vector, Candidate_Gene_log2Count, Candidate_TF_Peak_Binding, Candidate_Peak_Gene_looping,...
    T_A_prior, T_R_prior, T_prior_mean, T_prior_var, B_prior, B_prior_mean, B_prior_var, B_prob, L_prior, L_prior_mean, L_prior_var, L_prob, M, S, P, G, iteration_num)


iteration_seg=iteration_num/10;
A_sample=Candidate_Peak_log2Count;
R_sample=Candidate_Gene_log2Count;

B=B_prior;
L=L_prior;
T_A=T_A_prior;
T_R=T_R_prior;

B_state=Candidate_TF_Peak_Binding;
L_state=Candidate_Peak_Gene_looping;

B_state_frq=B_state;
L_state_frq=L_state;

alpha_A=1;
beta_A=1;
sigma_A_noise=var(A_sample-B_prior*T_prior_mean, 0, 'all');

alpha_R=1;
beta_R=1;
sigma_R_noise=var(R_sample-L_prior'*(B_prior*T_prior_mean), 0, 'all');


for i=1:iteration_num
    
    %********** scATAC-seq fitting and variable sampling  ****************
          
    %Step 1: TF activity sampling
    [T_A, T_R, T_sample] = TF_activity_T_sampling(scATAC_read_count_matrix, ATAC_cell_vector, A_sample, scRNA_read_count_matrix, RNA_cell_vector, R_sample, B, T_A, T_R, T_prior_mean, T_prior_var, sigma_A_noise, M, S, P, G);
    
    %Step 2: TF-peak binding weight sampling
    B = TF_peak_binding_B_sampling(scATAC_read_count_matrix, ATAC_cell_vector, A_sample, B, T_A, B_state, B_prior_mean, B_prior_var, sigma_A_noise, M, S, P, G);
    
    %Step 3: TF-peak binding state update
    [B_state, B] = TF_peak_binary_binding_B_state_sampling(scATAC_read_count_matrix, ATAC_cell_vector, A_sample, B, T_A, B_state, B_prior_mean, B_prior_var, B_prob, sigma_A_noise, M, S, P, G);
    
    %Step 4: ATAC fitting residue variance control    
    sigma_A_noise = 1/gamrnd(alpha_A+1/2,1/(beta_A+sum(sum((A_sample-B*T_sample).^2))/(2*P*S)));
    
       
    
    %********** scRNA-seq fitting and variable sampling  *****************      
    
    %Step 5: Peak-Gene looping weight sampling
    L = Peak_gene_looping_L_samping(scRNA_read_count_matrix, RNA_cell_vector, R_sample, L, B, T_R, L_state, L_prior_mean, L_prior_var, sigma_R_noise, M, S, P, G);
    
    %Step 6: Peak-Gene looping state update
    [L_state, L]=Peak_gene_binary_looping_L_state_samping(scRNA_read_count_matrix, RNA_cell_vector, R_sample, L, B, T_R, L_state, L_prior_mean, L_prior_var, L_prob, sigma_R_noise, M, S, P, G);
    
    %Step 7: RNA fitting residue variance control
    sigma_R_noise = 1/gamrnd(alpha_R+1/2,1/(beta_R+sum(sum((R_sample-L'*(B*T_sample)).^2))/(2*G*S)));
        
    
    
    %********** sample summary  *****************
        
    B_state_frq=B_state_frq+B_state;
    L_state_frq=L_state_frq+L_state;
    
    if mod(i,iteration_seg)==0
        fprintf('MAGICAL finished %d percent\n\n', 10*i/iteration_seg)
    end
    
       
end


Candidate_TF_Peak_Binding_prob=B_state_frq/iteration_num;
Candidate_Peak_Gene_Looping_prob=L_state_frq/iteration_num;
