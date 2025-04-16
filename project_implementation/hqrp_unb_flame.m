function [ A_out, t_out, s ] = hqrp_unb_flame( A )

	[ ~, n ] = size( A );
	t = zeros( n,1 );
	s = 1:n;

	[ ATL, ATR, ...
	  ABL, ABR ] = FLA_Part_2x2( A, ...
								 0, 0, 'FLA_TL' );
  
	[ tT, ...
	  tB ] = FLA_Part_2x1( t, ...
							0, 'FLA_TOP' );
	
	[ sL, sR ] = FLA_Part_1x2( s, ...
								0, 'FLA_LEFT' );
	while ( size( ATL, 2 ) < size( A, 2  ) )
		% Perform any column pivoting before partitioning
		[ s0, s1, s2 ] = FLA_Repart_1x2_to_1x3( sL, sR, 1, 'FLA_RIGHT' );
		pivot_idx = find_pivot(ABR);
		if pivot_idx ~= 1
			s1 = s2(pivot_idx - 1);
			s2(pivot_idx - 1) = s1;
			ABR(:, [1, pivot_idx]) = ABR(:, [pivot_idx, 1]);
		end
		[ sL, sR ] = FLA_Cont_with_1x3_to_1x2(s0, s1, s2, 'FLA_LEFT');

		[ A00,  a01,     A02,  ...
			a10t, alpha11, a12t, ...
			A20,  a21,     A22 ] = FLA_Repart_2x2_to_3x3( ATL, ATR, ...
														ABL, ABR, ...
														1, 1, 'FLA_BR' );
	
		[ t0, ...
			~, ...
			t2 ] = FLA_Repart_2x1_to_3x1( tT, ...
										tB, ...
										1, 'FLA_BOTTOM' );
	
		%------------------------------------------------------------%
	
		[ alpha11, ...
			a21, tau1 ] = housev( alpha11, ...
									a21 );
								
		w12t = ( a12t + a21' * A22 )/ tau1;
		
		a12t = a12t - w12t;
		A22  = A22 - a21 * w12t;
	
		%------------------------------------------------------------%
	
		[ ATL, ATR, ...
			ABL, ABR ] = FLA_Cont_with_3x3_to_2x2( A00,  a01,     A02,  ...
												a10t, alpha11, a12t, ...
												A20,  a21,     A22, ...
												'FLA_TL' );
	
		[ tT, ...
			tB ] = FLA_Cont_with_3x1_to_2x1( t0, ...
											tau1, ...
											t2, ...
											'FLA_TOP' );
  
	end
  
	A_out = [ ATL, ATR
			  ABL, ABR ];
		  
	t_out = [ tT, 
			  tB ];
  
  return
  
function [ idx ] = find_pivot(A)
	[~, idx ] = max(vecnorm(A, 2, 1));
return