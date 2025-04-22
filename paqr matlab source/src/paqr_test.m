% Test script for householder_poqr
clc;
clear;

% Parameters
m = 100;  % Number of rows
n = 50;  % Number of columns

% Generate random matrix A
A = randn(m, n);

% Choose orthogonalization routine
orth_qr2 = @orth_geqr2;  % or use orth_golub, orth_lapack, etc.

% Define is_deficient function (as in your test.m)
is_deficient = @(R, k, nrm, A, nrmA) abs(R(k,k)) < eps * nrmA;

% Choose norm function
nrm = @(A) norm(A, 2);
nrmA = nrm(A);

% Run POQR
[V, R, T, dead_cols] = householder_poqr(A, orth_qr2, is_deficient, nrm, nrmA);

% Reconstruct Q from V and T
Q = eye(m) - V * T * V';
Q = Q(:, ~dead_cols);


A_hat = Q * R;


% Compute error
error = norm(A_hat - A, 'fro') / norm(A, 'fro');

fprintf('Reconstruction error (Frobenius norm): %e\n', error);
