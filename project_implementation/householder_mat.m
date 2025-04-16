function [ H ] = householder_mat( u, tau )
    %{
        Input: u = [1, u2]'
        
        Output: H = I - 1/tau*u*u'
    %}
    n = length(u);
    H = eye(n) - 1/tau*(u*u');
end