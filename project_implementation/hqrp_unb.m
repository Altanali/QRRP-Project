function [ A, T, s ] = hqrp_unb(A, T, b)
%{
	Householder QR factorization with column pivoting
	[A, T, s] = HQR_UNB(A, T, b) returns the result of b iterations of 
	the Householder QR factorization with column pivoting of the matrix A.

	s contains the pivots used in the factorization.
%}
n = size(A, 2);
if b == -1
	b = n; % Default value for b
end
s = 1:n; % Initialize the pivot vector
for i = 1:b
	% Find the pivot column
	pivot = find_pivot(A(i: end, i: end));
	if pivot ~= i
		% Swap the columns in A and T
		A(:, [i, pivot]) = A(:, [pivot, i]);
		T(:, [i, pivot]) = T(:, [pivot, i]);
		s([i, pivot]) = s([pivot, i]); % Update the pivot vector
	end
	
	%Compute Householder Vector
	[A(i, i), A(i+1: end, i), T(i, i)] = housev(A(i, i), A(i+1: end, i));
end

function [ idx ] = find_pivot(A)
	[~, idx ] = max(vecnorm(A, 2, 1));
return