function [ A, T, s ] = hqrrp_blk( A, T, s, b, p )
	G = randn(b + p, size(A, 1));
	Y = G*A;

	[ ATL, ATR, ...
	  ABL, ABR ]  = FLA_Part_2x2( A, ...
	  							  0, 0, 'FLA_TL');

	[ TT, ...
	  TB ] = FLA_Part_2x1(T, 0, 'FLA_TOP');

	[ sL, sR ] = FLA_Part_1x2(s, 0, 'FLA_LEFT');

	[ YL, YR ] = FLA_Part_1x2(Y, 0, 'FLA_LEFT');
	[ GL, GR ] = FLA_Part_1x2(G, 0, 'FLA_LEFT');

	while (size(ATL, 2) < size(A, 2))
		b = min([b, size(ABR, 2)]);
		[ A00,  A01, A02,  ...
		  A10,  A11, A12, ...
		  A20,  A21, A22 ] = FLA_Repart_2x2_to_3x3( ATL, ATR, ...
														ABL, ABR, ...
														b, b, 'FLA_BR' );
		[ T0, ... 
		  T1, ..., 
		  T2 ] = FLA_Repart_2x1_to_3x1( TT, TB, b, 'FLA_BOTTOM' );
		
		% disp([sT]);
		% disp([sB]);
		[ s0, s1, s2 ] = FLA_Repart_1x2_to_1x3( sL, sR, b, 'FLA_RIGHT' );
		[ Y0, Y1, Y2 ] = FLA_Repart_1x2_to_1x3( YL, YR, b, 'FLA_RIGHT' );
		[ G0, G1, G2 ] = FLA_Repart_1x2_to_1x3( GL, GR, b, 'FLA_RIGHT' );

		s1 = DeterminePivots(Y1, Y2, b); %s1 is 1xb
		%Swap according to pivot selection
		[ A01, A02 ] = SwapCols(s1, A01, A02);
		[ A11, A12 ] = SwapCols(s1, A11, A12);
		[ A21, A22 ] = SwapCols(s1, A21, A22);
		
		%Perform HQRP with pivoted columns. 
		A_11_m = size(A11, 1);
		disp(["T1 here: ", size(T1)]);

		[A_out, T1, s1_prime] = hqrp_unb_flame([A11; A21], T1, b);
		disp(["T1 here: ", size(T1)]);

		A11 = A_out(1:A_11_m, :);
		A21 = A_out(A_11_m + 1:end, :);
 		[A01, ~] = SwapCols(s1_prime, A01, []);
		s1 = SwapCols(s1_prime, s1, []);

		%Update the rest of the matrix.
		U11 = tril(A11, -1) + eye(size(A11, 1));
		U21 = A21;
		R12 = A12;
		T11 = triu(T1(:, 1:b));
		W12 = T11' \ (A11' * A12 + A21' * A22);
		A12 = A12 - U11 * W12;
		A22 = A22 - U21 * W12;

		disp(["u11: ", size(U11)]);
		disp(["u21: ", size(U21)]);
		disp(["r12: ", size(R12)]);	
		disp(["t11: ", size(T11)]);
		disp(["w12: ", size(W12)]);
		disp(["a12: ", size(A12)]);
		disp(["a22: ", size(A22)]);
		disp(["y1: ", size(Y1)]);
		disp(["y2: ", size(Y2)]);
		disp(["g1: ", size(G1)]);
		disp(["g2: ", size(G2)]);
		%Update Y and G matrices.
		disp("------------------");

		[Y1, Y2] = SwapCols(s1, Y1, Y2);	
		[G1, G2] = SwapCols(s1, G1, G2);
		Y2 = Y2 - (G1 - ((G1 * U11 + G2 * U21) * inv(T11)) * U11') * R12;

		%Continue With...
		[ ATL, ATR, ... 
		  ABL, ABR ] = FLA_Cont_with_3x3_to_2x2( A00, A01, A02, ...
												A10, A11, A12, ...
												A20, A21, A22, ...
												'FLA_TR' );
		disp(["T0", size(T0)]);
		disp(["T1", size(T1)]);
		disp(["T2", size(T2)])

		[ TT, ... 
		  TB ] = FLA_Cont_with_3x1_to_2x1( T0, T1, T2, 'FLA_TOP' );

		[ sL, sR ] = FLA_Cont_with_1x3_to_1x2( s0, s1, s2, 'FLA_LEFT' );

		[ YL, YR ] = FLA_Cont_with_1x3_to_1x2( Y0, Y1, Y2, 'FLA_LEFT' );
		[ GL, GR ] = FLA_Cont_with_1x3_to_1x2( G0, G1, G2, 'FLA_LEFT' );
		
	end

	A = [ATL, ATR; 
	     ABL, ABR];
	T = [TT; TB];
	s = [sL, sR];
end


function [ s1 ] = DeterminePivots(Y1, Y2, b)
	[~, ~, s ] = hqrp_unb_flame(cat(2, Y1, Y2), -1, b);
	s1 = s( 1 : b );
end

function [ AL, AR ] = SwapCols( s1, AL, AR )
	% Permutes the columns of A01, A02, A11, A12, A21, A22 according to s1
	AL_n = size(AL, 2);
	A_temp = [AL, AR];
	perm = [s1, setdiff(1:size(A_temp, 2), s1)];
    A_temp = A_temp(:, perm);	
	AL = A_temp(:, 1 : AL_n);
	AR = A_temp(:, AL_n + 1 : end);
end