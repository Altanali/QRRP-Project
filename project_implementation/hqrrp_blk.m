function [ A, T, s ] = hqrrp_blk( A, T, s, b, p )
	if T == -1
		T = zeros(size(A, 2));
	end

	if s == -1
		s = 1:size(A, 2);
	end

	G = randn(b + p, size(A, 1));
	Y = G*A;

	[ ATL, ATR, ...
	  ABL, ABR ]  = FLA_Part_2x2( A, ...
	  							  0, 0, 'FLA_TL');

	% [ TT, ...
	%   TB ] = FLA_Part_2x1(T, 0, 'FLA_TOP');
	[ TTL, TTR, ... 
	  TBL, TBR ] = FLA_Part_2x2( T, ... 
	  							  0, 0, 'FLA_TL' );
	

	[ sL, sR ] = FLA_Part_1x2(s, 0, 'FLA_LEFT');
	[ YL, YR ] = FLA_Part_1x2(Y, 0, 'FLA_LEFT');
	[ GL, GR ] = FLA_Part_1x2(G, 0, 'FLA_LEFT');

	while (size(ATL, 2) < size(A, 2) && size(ATL, 1) < size(A, 1))
		b = min([b, size(ABR, 2)]);
		[ A00,  A01, A02,  ...
		  A10,  A11, A12, ...
		  A20,  A21, A22 ] = FLA_Repart_2x2_to_3x3( ATL, ATR, ...
														ABL, ABR, ...
														b, b, 'FLA_BR' );
		% [ T0, ... 
		%   T1, ..., 
		%   T2 ] = FLA_Repart_2x1_to_3x1( TT, TB, b, 'FLA_BOTTOM' );
		[ T00,  T01, T02,  ...
		  T10,  T11, T12, ...
		  T20,  T21, T22 ] = FLA_Repart_2x2_to_3x3( TTL, TTR, ... 
														TBL, TBR, ... 
														b, b, 'FLA_BR' );
		[ s0, s1, s2 ] = FLA_Repart_1x2_to_1x3( sL, sR, b, 'FLA_RIGHT' );
		[ Y0, Y1, Y2 ] = FLA_Repart_1x2_to_1x3( YL, YR, b, 'FLA_RIGHT' );
		[ G0, G1, G2 ] = FLA_Repart_1x2_to_1x3( GL, GR, b, 'FLA_RIGHT' );
        [~, ~, s_pivots] = hqrp_unb_flame(YR, -1, -1, b);



		%Swap according to pivot selection
		[ A01, A02 ] = SwapCols(s_pivots, A01, A02);
		[ A11, A12 ] = SwapCols(s_pivots, A11, A12);
		[ A21, A22 ] = SwapCols(s_pivots, A21, A22);
        [ Y1, Y2 ] = SwapCols(s_pivots, Y1, Y2);
		[s1, s2] = SwapCols(s_pivots, s1, s2);

		%Perform HQRP with pivoted columns. 
		[A_out, T11, s1_prime] = hqrp_unb_flame([A11; A21], T11, -1, -1);
        
		A11 = A_out(1:b, :);
		A21 = A_out(b + 1:end, :);

        %Apply permutations to A01 and Y1
 		[A01, ~] = SwapCols(s1_prime, A01, []);
        [Y1, ~] = SwapCols(s1_prime, Y1, []);
		[s1, ~] = SwapCols(s1_prime, s1, []);
        

		%Update the rest of the matrix.
		U11 = tril(A11, -1) + eye(size(A11, 1));
		U21 = A21;
		R12 = A12;
		T11_T = triu(T11(:, 1:b));
		W12 = inv(T11_T') * (U11' * A12 + U21' * A22);
		A12 = A12 - U11 * W12;
		A22 = A22 - U21 * W12;

		%Update Y and G matrices.
		if size(Y2, 2) > 0
			Y2 = Y2 - (G1 - ((G1 * U11 + G2 * U21) * inv(T11_T)) * U11') * R12;
		end

		%Continue With...
		[ ATL, ATR, ... 
		  ABL, ABR ] = FLA_Cont_with_3x3_to_2x2( A00, A01, A02, ...
												A10, A11, A12, ...
												A20, A21, A22, ...
												'FLA_TL' );

		% [ TT, ... 
		%   TB ] = FLA_Cont_with_3x1_to_2x1( T0, T1, T2, 'FLA_TOP' );
		[ TTL, TTR, ... 
		  TBL, TBR ] = FLA_Cont_with_3x3_to_2x2( T00, T01, T02, ...
												T10, T11, T12, ...
												T20, T21, T22, ...
												'FLA_TL' );
		[ sL, sR ] = FLA_Cont_with_1x3_to_1x2( s0, s1, s2, 'FLA_LEFT' );
		[ YL, YR ] = FLA_Cont_with_1x3_to_1x2( Y0, Y1, Y2, 'FLA_LEFT' );
		[ GL, GR ] = FLA_Cont_with_1x3_to_1x2( G0, G1, G2, 'FLA_LEFT' );
		
	end

	A = [ATL, ATR; 
	     ABL, ABR];
	% T = [TT; TB];
	T = [TTL, TTR; 
	     TBL, TBR];
	s = [sL, sR];
end



function [ AL, AR ] = SwapCols( perm, AL, AR )
    % len(perm) = col(AL) + col(AR)
	% Permutes the columns of A01, A02, A11, A12, A21, A22 according to s1
	AL_n = size(AL, 2);
	A_temp = [AL, AR];
    A_temp = A_temp(:, perm);	
	AL = A_temp(:, 1 : AL_n);
	AR = A_temp(:, AL_n + 1 : end);
end