function [B_state_new, B_new] = TF_peak_binary_binding_B_state_sampling(A, B, T, B_state, B_mean, B_var, B_prob, sigma_A_noise)


F=size(B,1);
M=size(B,2);
S=size(A,1);
TF_index=randperm(M);
for f=1:F
    temp=B(f,:)*T;
    for i=1:M
        m=TF_index(i);
        temp_var=(T(m,:)*T(m,:)'*B_var(m)/S+sigma_A_noise);
        
        if B_state(f,m)==1
            
            mean_B=((A(f,:)-temp+B(f,m)*T(m,:))*T(m,:)'*B_var(m)/S+B_mean(f,m)*sigma_A_noise)/temp_var;
            vairance_B=B_var(m)*sigma_A_noise/temp_var;
            
            post_b1 = exp(-(B(f,m)*1-mean_B)^2/(2*vairance_B))*(B_prob(f,m)+0.25);
            post_b0 = exp(-(B(f,m)*0-mean_B)^2/(2*vairance_B))*(1-B_prob(f,m)+0.25);
            
            P1=post_b1/(post_b1+post_b0);
            threshold_c=rand;
            
            if P1>=threshold_c
                %hold the value
                B(f,m) = B(f,m);
                B_state(f,m) = B_state(f,m);
            else
                %flip 1 to 0
                B(f,m) = 0;
                B_state(f,m) = 0;
            end
            
        elseif B_state(f,m)==0 && B_mean(f,m)~=0
            
            mean_B=((A(f,:)-temp)*T(m,:)'*B_var(m)/S+B_mean(f,m)*sigma_A_noise)/temp_var;
            vairance_B=B_var(m)*sigma_A_noise/temp_var;
            
            B_temp=randn*sqrt(vairance_B)+mean_B;
            
            post_b1 = exp(-(B_temp*1-mean_B)^2/(2*vairance_B))*(B_prob(f,m)+0.25);
            post_b0 = exp(-(B_temp*0-mean_B)^2/(2*vairance_B))*(1-B_prob(f,m)+0.25);
            
            P1=post_b1/(post_b1+post_b0);
            threshold_c=rand;
            
            if P1>=threshold_c
                %flip 0 to 1 and assign the new value
                B(f,m) = B_temp;
                B_state(f,m) = 1;
            else
                %hold the value
                B(f,m) = B(f,m);
                B_state(f,m) = B_state(f,m);
            end
            
        else
            %no motif mapping positions are always 0
            B(f,m) = 0;
            B_state(f,m) = 0;
            
        end
    end
end

B_new=B;
B_state_new=B_state;
