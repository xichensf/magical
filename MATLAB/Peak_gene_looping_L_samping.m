function L_new = Peak_gene_looping_L_samping(R, L, A_estimate, L_state, L_mean, L_var, sigma_R_noise)
% L_new=L;
G=size(R,1);
F=size(L,1);
S=size(R,2);
F_index=randperm(F);
for i=1:F
    f=F_index(i);
    temp_var=A_estimate(f,:)*A_estimate(f,:)'*L_var/S+sigma_R_noise;%1 by 1 number
    mean_L=(A_estimate(f,:)*(R-L'*A_estimate+L(f,:)'*A_estimate(f,:))'*L_var/S+L_mean(f,:)*sigma_R_noise)/temp_var;%1 by G vector
    vairance_L=L_var*sigma_R_noise/temp_var;%1 by 1 number
    L(f,:)=(randn(1,G)*sqrt(vairance_L)+mean_L).*L_state(f,:);%1 by G vector
end
L_new=L;
