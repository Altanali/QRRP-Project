results = struct('test_paqr', [], 'test_hqrrp', [], 'test_pa_hqrrp', [], 'test_hqr', [], 'test_hqrp', []);
m = 1000
for rk = 100:100:m
     % Generate random matrices and vectors
    A = zeros(m, m); 
    %Create a rank rk matrix via a sum of rk outer products. 
    for i = 1:rk
        x_temp = rand(m, 1);
        A = A + x_temp * x_temp';
    end
    x = rand(m, 1);
    % Run each test function and store the results
    [time, error] = test_paqr(A, x);
    results.test_paqr = [results.test_paqr; rk, time, error];
    
    [time, error] = test_hqrrp(A, x);
    results.test_hqrrp = [results.test_hqrrp; rk, time, error];
    
    [time, error] = test_pa_hqrrp(A, x);
    results.test_pa_hqrrp = [results.test_pa_hqrrp; rk, time, error];
    
    [time, error] = test_hqr(A, x);
    results.test_hqr = [results.test_hqr; rk, time, error];
    
    [time, error] = test_hqrp(A, x);
    results.test_hqrp = [results.test_hqrp; rk, time, error];
end

% Save the results to a .mat file for later use
save('project_implementation/results.mat', 'results');
%%

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
    Q = Q(:, ~dead_cols);
    r = size(R, 2);
    R = R(1:r, :);
    x_partial = R \ (Q' * b);
    x_hat = zeros(n, 1);
    x_hat(~dead_cols) = x_partial;
    
    
    error = norm(A' * (A*x_hat - b), 2) / (norm(A, 2)^2);
    disp(["paqr error: ", error]);
end
%%
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
    disp(["hqrrp_blk error:  ", error]);
end
%%
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
    disp(["pa_hqrrp_blk error:  ", error]);
end
%%
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
    disp(["hqr error: ", error]);
end
%%
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
    disp(["hqrp error: ", error]);
end