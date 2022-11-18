function [B_state, B] = TF_peak_binary_binding_B_state_sampling(scATAC_read_count_matrix, ATAC_cell_vector, A_sample, B, T_A,...
    B_state, B_mean, B_var, B_prob, sigma_A_noise, M, S, P, G)
  
T_sample=zeros(M,S);
for s=1:S
    T_sample(:,s)=mean(T_A(:,ATAC_cell_vector==s), 2);
end
    
  
TF_index=randperm(M);
for p=1:P
    temp=B(p,:)*T_sample;
    for i=1:M
        m=TF_index(i);
        temp_var=(T_sample(m,:)*T_sample(m,:)'*B_var(m)/S+sigma_A_noise);
        
        if B_state(p,m)==1
            
            mean_B=((A_sample(p,:)-temp+B(p,m)*T_sample(m,:))*T_sample(m,:)'*B_var(m)/S+B_mean(p,m)*sigma_A_noise)/temp_var;
            vairance_B=B_var(m)*sigma_A_noise/temp_var;
            
            post_b1 = exp(-(B(p,m)*1-mean_B)^2/(2*vairance_B))*(B_prob(p,m)+0.25)+1e-6;
            post_b0 = exp(-(B(p,m)*0-mean_B)^2/(2*vairance_B))*(1-B_prob(p,m)+0.25)+1e-6;
            
            P1=post_b1/(post_b1+post_b0);
            threshold_c=rand;
            
            if P1>=threshold_c
                %hold the value
                B(p,m) = B(p,m);
                B_state(p,m) = B_state(p,m);
            else
                %flip 1 to 0
                B(p,m) = 0;
                B_state(p,m) = 0;
            end
            
        elseif B_state(p,m)==0 && B_mean(p,m)~=0
            
            mean_B=((A_sample(p,:)-temp)*T_sample(m,:)'*B_var(m)/S+B_mean(p,m)*sigma_A_noise)/temp_var;
            vairance_B=B_var(m)*sigma_A_noise/temp_var;
            bb=randn;
            
            if bb>3
                bb=3;
            end
            
            if bb<-3
                bb=-3;
            end
            
            B_temp=bb*sqrt(vairance_B)+mean_B;
            
            post_b1 = exp(-(B_temp*1-mean_B)^2/(2*vairance_B))*(B_prob(p,m)+0.25)+1e-6;
            post_b0 = exp(-(B_temp*0-mean_B)^2/(2*vairance_B))*(1-B_prob(p,m)+0.25)+1e-6;
            
            P1=post_b1/(post_b1+post_b0);
            threshold_c=rand;
            
            if P1>=threshold_c
                %flip 0 to 1 and assign the new value
                B(p,m) = B_temp;
                B_state(p,m) = 1;
            else
                %hold the value
                B(p,m) = B(p,m);
                B_state(p,m) = B_state(p,m);
            end
            
        else
            %no motif mapping positions are always 0
            B(p,m) = 0;
            B_state(p,m) = 0;
            
        end
    end
end
