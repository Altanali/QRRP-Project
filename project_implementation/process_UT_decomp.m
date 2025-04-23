function [A_perm, Q, R] = process_UT_decomp(A_orig, A_out, T_out, perm)
    [m, n] = size(A_orig);
    k = min([m, n]);
    A_perm = A_orig(:, perm);
    R = triu(A_out);
    U = eye(m, n) + tril(A_out, -1);
    U = U(:, 1:k);
    Q = eye(m, m) - U*inv(T_out)*U';

end

