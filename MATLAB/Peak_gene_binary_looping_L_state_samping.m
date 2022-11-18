function  [L_state, L]=Peak_gene_binary_looping_L_state_samping(scRNA_read_count_matrix, RNA_cell_vector, R_sample, L, B, T_R, L_state, L_mean, L_var, L_prob, sigma_R_noise, M, S, P, G)

T_sample=zeros(M,S);

for s=1:S
    T_sample(:,s)=mean(T_R(:,RNA_cell_vector==s), 2);
end
 

A_estimate=B*T_sample;

for g=1:G
    temp=L(:,g)'*A_estimate;
    for p=1:P
        temp_var=A_estimate(p,:)*A_estimate(p,:)'*L_var/S+sigma_R_noise;
        
        if L_state(p,g)==1
            
            mean_L=((R_sample(g,:)-temp+L(p,g)*A_estimate(p,:))*A_estimate(p,:)'*L_var/S+L_mean(p,g)*sigma_R_noise)/temp_var;
            vairance_L=L_var*sigma_R_noise/temp_var;
            
            post_l1 = exp(-(L(p,g)*1-mean_L)^2/(2*vairance_L))*(L_prob(p,g)+0.1)+1e-6;
            post_l0 = exp(-(L(p,g)*0-mean_L)^2/(2*vairance_L))*(1-L_prob(p,g)+0.1)+1e-6;
            
            P1=post_l1/(post_l1+post_l0);
            threshold_l=rand;
            
            if P1>=threshold_l
                %hold the value
                L(p,g) = L(p,g);
                L_state(p,g) = L_state(p,g);
            else
                %flip 1 to 0
                L(p,g) = 0;
                L_state(p,g) = 0;
            end
            
            
        elseif L_state(p,g)==0 && L_mean(p,g)~=0
            
            mean_L=((R_sample(g,:)-temp)*A_estimate(p,:)'*L_var/S+L_mean(p,g)*sigma_R_noise)/temp_var;
            vairance_L=L_var*sigma_R_noise/temp_var;
            
            ll=randn;
            
            if ll>3
                ll=3;
            end
            
            if ll<-3
                ll=-3;
            end
            
            L_temp=ll*sqrt(vairance_L)+mean_L;%propose a new weight
            
            post_l1 = exp(-(L_temp*1-mean_L)^2/(2*vairance_L))*(L_prob(p,g)+0.1)+1e-6;
            post_l0 = exp(-(L_temp*0-mean_L)^2/(2*vairance_L))*(1-L_prob(p,g)+0.1)+1e-6;
            
            P1=post_l1/(post_l1+post_l0);
            threshold_l=rand;
            
            if P1>=threshold_l
                %flip 0 to 1 and assign a new value
                L(p,g) = L_temp;
                L_state(p,g) = 1;
            else
                %hold the value 0
                L(p,g) = L(p,g);
                L_state(p,g) = L_state(p,g);
            end
            
        else
            %peaks and genes more than 500K are not evaluated
            L(p,g) = 0;
            L_state(p,g) = 0;
            
        end
    end
end