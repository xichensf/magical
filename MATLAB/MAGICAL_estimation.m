function [Candidate_TF_Peak_Binding_prob, Candidate_Peak_Gene_Looping_prob]=...
    MAGICAL_estimation(A, R, B_state, L_state, T_prior, T_mean, T_var, B_prior, B_mean, B_var, B_prob, L_prior, L_mean, L_var, L_prob, S, P, G, iteration_num)


iteration_seg=iteration_num/10;


T=T_prior;
B=B_prior;
L=L_prior;


B_state_frq=B_state;
L_state_frq=L_state;

alpha_A=1;
beta_A=1;
ATAC_fitting_residue=A-B*T;
sigma_A_noise=var(ATAC_fitting_residue(:));

% RNA_fitting_residue=R-L'*(B*P);
% sigma_R_noise=var(RNA_fitting_residue(:));
alpha_R=1;
beta_R=1;
RNA_fitting_residue=R-L'*A;
sigma_R_noise=var(RNA_fitting_residue(:));


for i=1:iteration_num
    
    %********** scATAC-seq fitting and variable sampling  ****************
    
    %Step 1: TF activity sampling
    T = TF_activity_T_sampling(A, B, T, T_mean, T_var, sigma_A_noise);
    
    %Step 2: TF-peak binding weight sampling
    B = TF_peak_binding_B_sampling(A, B, T, B_state, B_mean, B_var, sigma_A_noise);
    
    %Step 3: TF-peak binding state update
    [B_state, B] = TF_peak_binary_binding_B_state_sampling(A, B, T, B_state, B_mean, B_var, B_prob, sigma_A_noise);
    
    %Step 4: ATAC fitting residue variance control
    ATAC_fitting_residue=A-B*T;
    
    sigma_A_noise = 1/gamrnd(alpha_A+1/2,1/(beta_A+sum(sum(ATAC_fitting_residue.^2))/(2*P*S)));
    % sigma_A_noise = (beta_A+sum(sum(ATAC_fitting_residue.^2))/(2*F*S))/chi2rnd(2*alpha_A+1);
    noise_A_vector(i)=sigma_A_noise;
    
    %********** scRNA-seq fitting and variable sampling  *****************
    A_estimate=B*T;% or true A
    
    %Step 5: Peak-Gene looping weight sampling
    L = Peak_gene_looping_L_samping(R, L, A_estimate, L_state, L_mean, L_var, sigma_R_noise);
    
    %Step 6: Peak-Gene looping state update
    [L_state, L]=Peak_gene_binary_looping_L_state_samping(R, L, A_estimate, L_state, L_mean, L_var, L_prob, sigma_R_noise);
    
    %Step 7: RNA fitting residue variance control
    RNA_fitting_residue=R-L'*A_estimate;
    sigma_R_noise = 1/gamrnd(alpha_R+1/2,1/(beta_R+sum(sum(RNA_fitting_residue.^2))/(2*G*S)));
    %         sigma_R_noise = (beta_R+sum(sum(RNA_fitting_residue.^2))/(2*G*S))/chi2rnd(2*alpha_R+1);
    noise_B_vector(i) = sigma_R_noise;
    
    %Step 8: Sample Summary
    B_state_frq=B_state_frq+B_state;
    L_state_frq=L_state_frq+L_state;
    
    if mod(i,iteration_seg)==0
        fprintf('MAGICAL finished %d percent\n\n', 10*i/iteration_seg)
    end
end


Candidate_TF_Peak_Binding_prob=B_state_frq/iteration_num;
Candidate_Peak_Gene_Looping_prob=L_state_frq/iteration_num;
