function [L_state_new, L_new]=Peak_gene_binary_looping_L_state_samping(R, L, A_estimate, L_state, L_mean, L_var, L_prob, sigma_R_noise)
% L_new=L;
% L_state_new=L_state;

G=size(L,2);
F=size(L,1);
S=size(R,2);

for g=1:G
    temp=L(:,g)'*A_estimate;

    for f=1:F
        temp_var=A_estimate(f,:)*A_estimate(f,:)'*L_var/S+sigma_R_noise;
        
        if L_state(f,g)==1
            
            mean_L=((R(g,:)-temp+L(f,g)*A_estimate(f,:))*A_estimate(f,:)'*L_var/S+L_mean(f,g)*sigma_R_noise)/temp_var;
            vairance_L=L_var*sigma_R_noise/temp_var;
            
            post_l1 = exp(-(L(f,g)*1-mean_L)^2/(2*vairance_L))*(L_prob(f,g)+0.1);
            post_l0 = exp(-(L(f,g)*0-mean_L)^2/(2*vairance_L))*(1-L_prob(f,g)+0.1);
            
            P1=post_l1/(post_l1+post_l0);
            threshold_l=rand;
            
            if P1>=threshold_l
                %hold the value
                L(f,g) = L(f,g);
                L_state(f,g) = L_state(f,g);
            else
                %flip 1 to 0
                L(f,g) = 0;
                L_state(f,g) = 0;
            end
            
            
        elseif L_state(f,g)==0 && L_mean(f,g)~=0
            
            mean_L=((R(g,:)-temp)*A_estimate(f,:)'*L_var/S+L_mean(f,g)*sigma_R_noise)/temp_var;
            vairance_L=L_var*sigma_R_noise/temp_var;
            
            L_temp=randn*sqrt(vairance_L)+mean_L;%propose a new weight
            
            post_l1 = exp(-(L_temp*1-mean_L)^2/(2*vairance_L))*(L_prob(f,g)+0.1);
            post_l0 = exp(-(L_temp*0-mean_L)^2/(2*vairance_L))*(1-L_prob(f,g)+0.1);
            
            P1=post_l1/(post_l1+post_l0);
            threshold_l=rand;
            
            if P1>=threshold_l
                %flip 0 to 1 and assign a new value
                L(f,g) = L_temp;
                L_state(f,g) = 1;
            else
                %hold the value 0
                L(f,g) = L(f,g);
                L_state(f,g) = L_state(f,g);
            end
            
        else
            %peaks and genes more than 500K are not evaluated
            L(f,g) = 0;
            L_state(f,g) = 0;
            
        end
    end
end

L_new=L;
L_state_new=L_state;
