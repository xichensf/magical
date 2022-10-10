function [B_state_new, B_new] = TF_peak_binary_binding_B_state_sampling(A, B, P, B_state, B_mean, B_var, B_prob, sigma_A_noise)

% B_state_new=B_state;
% B_new=B_state;
F=size(B,1);
T=size(B,2);
S=size(A,1);
TF_index=randperm(T);
for f=1:F
    temp=B(f,:)*P;
    for i=1:T
        t=TF_index(i);
        temp_var=(P(t,:)*P(t,:)'*B_var(t)/S+sigma_A_noise);
        
        if B_state(f,t)==1
            
            mean_B=((A(f,:)-temp+B(f,t)*P(t,:))*P(t,:)'*B_var(t)/S+B_mean(f,t)*sigma_A_noise)/temp_var;
            vairance_B=B_var(t)*sigma_A_noise/temp_var;
            
            post_b1 = exp(-(B(f,t)*1-mean_B)^2/(2*vairance_B))*(B_prob(f,t)+0.25);
            post_b0 = exp(-(B(f,t)*0-mean_B)^2/(2*vairance_B))*(1-B_prob(f,t)+0.25);
            
            P1=post_b1/(post_b1+post_b0);
            threshold_c=rand;
            
            if P1>=threshold_c
                %hold the value
                B(f,t) = B(f,t);
                B_state(f,t) = B_state(f,t);
            else
                %flip 1 to 0
                B(f,t) = 0;
                B_state(f,t) = 0;
            end
            
        elseif B_state(f,t)==0 && B_mean(f,t)~=0
            
            mean_B=((A(f,:)-temp)*P(t,:)'*B_var(t)/S+B_mean(f,t)*sigma_A_noise)/temp_var;
            vairance_B=B_var(t)*sigma_A_noise/temp_var;
            
            B_temp=randn*sqrt(vairance_B)+mean_B;
            
            post_b1 = exp(-(B_temp*1-mean_B)^2/(2*vairance_B))*(B_prob(f,t)+0.25);
            post_b0 = exp(-(B_temp*0-mean_B)^2/(2*vairance_B))*(1-B_prob(f,t)+0.25);
            
            P1=post_b1/(post_b1+post_b0);
            threshold_c=rand;
            
            if P1>=threshold_c
                %flip 0 to 1 and assign the new value
                B(f,t) = B_temp;
                B_state(f,t) = 1;
            else
                %hold the value
                B(f,t) = B(f,t);
                B_state(f,t) = B_state(f,t);
            end
            
        else
            %no motif mapping positions are always 0
            B(f,t) = 0;
            B_state(f,t) = 0;
            
        end
    end
end

B_state_new=B_state;
B_new=B_state;