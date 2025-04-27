function [error, time ] = test_qr_algo(A, x, algo_name)
    if strcmp(algo_name, "paqr")
        [error, time] = test_paqr(A, x);
    elseif strcmp(algo_name, "hqrrp")
        [error, time] = test_hqrrp(A, x);
    elseif strcmp(algo_name, "pa_hqrrp")
        [error, time] = test_pa_hqrrp(A, x);
    elseif strcmp(algo_name, "hqr")
        [error, time] = test_hqr(A, x);
    elseif strcmp(algo_name, "hqrp")
        [error, time] = test_hqrp(A, x);
    end
end

function [error, time] = test_paqr(A, x) 
    [m, n] = size(A);
    b = A * x;
    % Choose orthogonalization routine
    orth_qr2 = @orth_geqr2;  % or use orth_golub, orth_lapack, etc.
    
    % Define is_deficient function (as in your test.m)
    is_deficient = @(R, k, nrm, A, nrmA) abs(R(k,k)) < eps * nrmA;
    
    % Choose norm function
    nrm = @(A) norm(A, 2);
    nrmA = nrm(A);
    % Run POQR
    tic;
    [V, R, T, dead_cols] = householder_poqr(A, orth_qr2, is_deficient, nrm, nrmA);
    time = toc;
    Q = eye(m) - V * T * V';
    x_partial = R \ (Q' * b);
    x_hat = zeros(n, 1);
    x_hat(~dead_cols) = x_partial;
    
    
    error = norm(A' * (A*x_hat - b), 2) / (norm(A, 2)^2);
    %disp(["paqr error: ", error]);
end

function [error, time] = test_hqrrp(A, x)
    tic
    [A_out, T_out, s] = hqrrp_blk(A, -1, -1, 5, 10);
    time = toc;
    tau = diag(T_out);
    [A_perm, Q1, R1] = process_VT_decomp(A, A_out, tau, s);
    b_hat = A_perm * x;
    Q = Q1;
    x_hat = R1 \ (Q' * b_hat);
    error = norm(A_perm' * (A_perm*x_hat - b_hat), 2) / (norm(A_perm, 2)^2);
    %disp(["hqrrp_blk error:  ", error]);
end

function [error, time] = test_pa_hqrrp(A, x)
    %{
        Test Pivot Avoiding Householder QR Randomizaton for Column
        Pivoting.
        Error calculation: Calculate APx, where P is the matrix 
        denoting the permutation and and deletion of columns of A
        according to to pa_hqrrp_blk.
        Move all dropped columns to the end of the permuted matrix in 
        an arbitrary order. 
    %}
    block_size = 64;
    oversample_size = 10;
    tic
    [A_out, T_out, s, dead_cols] = pa_hqrrp_blk(A, -1, -1, block_size, oversample_size, 1);
    time = toc;
    num_columns_dropped = sum(dead_cols);
    tau = diag(T_out);
    perm = cat(2, s, setdiff(1:size(A, 2), s));
    [A_perm, Q, R] = process_VT_decomp(A, A_out, tau, perm);
    R = R(1:size(s, 2), :);
    
    b_hat = A_perm * x;
    x_partial = R \ (Q' * b_hat);
    x_hat = cat(1, x_partial, zeros(num_columns_dropped, 1));
    error = norm(A_perm' * (A_perm*x_hat - b_hat), 2) / (norm(A_perm, 2)^2);
    %disp(["pa_hqrrp_blk error:  ", error]);
end

function [error, time] = test_hqr(A, x)
    %{
        Test the unblocked householder QR transformation without column
        pivoting.
    %}
    tic 
    [A_out, T_out, s] = hqrp_unb_flame(A, -1, -1, -1, 0);
    time = toc;
    [~, Q, R] = process_UT_decomp(A, A_out, T_out, s);
    b = A * x;
    x_hat = R \ (Q' * b);
    error = norm(A' * (A*x_hat - b), 2) / (norm(A, 2)^2);
    %disp(["hqr error: ", error]);
end

function [error, time] = test_hqrp(A, x)
    %{
        Test the unblocked householder QR transformation with column
        pivoting.
    %}
    tic 
    [A_out, T_out, s] = hqrp_unb_flame(A, -1, -1, -1, 1);
    time = toc;
    [A_perm, Q, R] = process_UT_decomp(A, A_out, T_out, s);
    b = A * x;
    x_hat = R \ (Q' * b);
    error = norm(A_perm' * (A_perm*x_hat - b), 2) / (norm(A_perm, 2)^2);
    %disp(["hqrp error: ", error]);
end