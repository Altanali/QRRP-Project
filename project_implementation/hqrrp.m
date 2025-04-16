function [A, T, s] = hqrrp(A, T, s, b, p)
    G = randn(b + p, size(A, 2));

    Y = G * A;

    % Partition
    A_TL = []; A_TR = []; A_BL = []; A_BR = A
    T_T = []; T_B = T;
    s_T = []; s_B = s;
    Y_L = []; Y_R = Y;
    while size(A_TL, 1) < size(A, 1)
        % Determine block size
        b = min(b, size(A_BR, 2));

        % Repartition
        A_00 = A_TL; A_01 = A_TR(:, 1:b); A_02 = A_TR(:, b + 1:end);
        A_10 = A_BL(1:b, :); A_11 = A_BR(1:b, 1:b); A_12 = A_BR(1:b, b + 1:end);
        A_20 = A_BL(b + 1:end, :); A_21 = A_BR(b + 1:end, 1:b); A_22 = A_BR(b + 1:end, b + 1:end);
        T_0 = T_T; T_1 = T_B(1:b, :); T_2 = T_B(b + 1:end, :);
        s_0 = s_T; s_1 = s_B(1:b); s_2 = s_B(b + 1:end);
        Y_0 = Y_L; Y_1 = Y_R(:, 1:b); Y_2 = Y_R(:, b + 1:end);

        s_1 = DeterminePivots([Y_1, Y_2], b);

        A_blk = [A_01, A_02];
        perm = [s_1, setdiff(1:size(A_blk, 2), s_1)];
        A_blk = A_blk(:, perm);
        A_01 = A_blk(:, 1:length(s_1));
        A_02 = A_blk(:, length(s_1)+1:end);
        A_blk = [A_11, A_12];
        perm = [s_1, setdiff(1:size(A_blk, 2), s_1)];
        A_blk = A_blk(:, perm);
        A_11 = A_blk(:, 1:length(s_1));
        A_12 = A_blk(:, length(s_1)+1:end);
        A_blk = [A_21, A_22];
        perm = [s_1, setdiff(1:size(A_blk, 2), s_1)];
        A_blk = A_blk(:, perm);
        A_21 = A_blk(:, 1:length(s_1));
        A_22 = A_blk(:, length(s_1)+1:end);

        [A_blk, T_1, s_1prime] = hqrp_unb([A_11; A_21], T_1, b);
        A_11 = A_blk(1:b, :);
        A_21 = A_blk(b + 1:end, :);

        A_01 = A_01(:, [s_1prime, setdiff(1:size(A_01, 2), s_1prime)]);

        s_1 = s_1prime;

        U_11 = tril(A_11, -1) + eye(size(A_11, 1));
        U_21 = A_21
        W_12 = inv(T_1') * (U_11' * A_12 + U_21' * A_22);

        A_12 = A_12 - U_11 * W_12;
        A_22 = A_22 - U_21 * W_12;

        Y_blk = [Y_1, Y_2];
        perm = [s_1, setdiff(1:size(Y_blk, 2), s_1)];
        Y_blk = Y_blk(:, perm);
        Y_1 = Y_blk(:, 1:length(s_1));
        Y_2 = Y_blk(:, length(s_1)+1:end);

        R_12 = A_12
        Y_2 = Y_2 - ()

    end
end

function s = DeterminePivots(Y, b)
    [~, ~, p] = qr(Y, 'vector');   % p is a permutation of 1:n
    s = p(1:b);                    % select top-b pivot indices
end