% Test script for householder_poqr
clc;
clear;

% Parameters
m = 100;  % Number of rows
n = 50;  % Number of columns

% Generate random matrices and vectors
A = randn(m, 20) * randn(20, n); % Random matrix A
x_hat = randn(n, 1);
b = A * x_hat;


% Choose orthogonalization routine
orth_qr2 = @orth_geqr2;  % or use orth_golub, orth_lapack, etc.

% Define is_deficient function (as in your test.m)
is_deficient = @(R, k, nrm, A, nrmA) abs(R(k,k)) < eps * nrmA;

% Choose norm function
nrm = @(A) norm(A, 2);
nrmA = nrm(A);

% Run POQR
[V1, R1, T1, dead_cols] = householder_poqr(A, orth_qr2, is_deficient, nrm, nrmA);

Q1 = eye(m) - V1 * T1 * V1';
Q1 = Q1(:, ~dead_cols);
r = size(R1, 2);

R1 = R1(1:r, :);


[V2, R2, T2] = householder_qr(R1, orth_qr2);
Q2 = eye(r) - V2 * T2 * V2';

Q = Q1 * Q2;
x_partial = R2 \ (Q' * b);

x = zeros(n, 1);
x(~dead_cols) = x_partial;


error = norm(A' * (A*x - b), 2) / (norm(A, 2)^2);
disp(error);
