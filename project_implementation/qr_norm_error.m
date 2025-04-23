function [error] = qr_norm_error(A_orig, Q, R, perm, nrm)
    if isempty(nrm)
        nrm = @(A) norm(A, "fro");
    end
    
    if isempty(perm)
        perm = 1:size(A_orig, 2);
    end

    A_hat = Q*R;
    A_perm = A_orig(:, perm);
    error = nrm((A_perm - A_hat))/nrm(A_perm);
end

