function [A_perm, Q, R] = process_VT_decomp(A_orig, A_out, tau, perm)
    Q = FormQ(A_out, tau);
    R = triu(A_out);
    A_perm = A_orig(:, perm);
end

