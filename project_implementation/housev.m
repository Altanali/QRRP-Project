function [ rho, u2, tau ] = housev(chi1, x2)
    %{
        Input: x = [chi1, x2]'
        
        Output: 
            - u = [1, u2]': reflector
            - tau: normalizing factor for u
            - rho (diagonal entry) 
            
    %}
    chi2 = norm (x2);
    alpha = norm ([chi1, chi2]); % |x|_2
    rho = -sign(chi1)*alpha;
    v1 = chi1 - rho;
    u2 = x2/v1;
    chi2 = x2/abs(v1); % |u_2|_2
    tau = (1 + chi2^2)/2;
return