function [ A_out, T_out, s ] = hqrp_unb_flame( A , T , s, num_iter, pivot_true)
	[ ~, n ] = size( A );

	if num_iter == -1
		num_iter = (size(A, 2));
	end 
    
    if s == -1
        s = 1:n;
    end

	if T == -1
		T = zeros(size(A, 2));
	end 

	
	% t = zeros( n,1 );
	[ ATL, ATR, ...
	  ABL, ABR ] = FLA_Part_2x2( A, ...
								 0, 0, 'FLA_TL' );
	[ TTL, TTR, ..., 
	  TBL, TBR ] = FLA_Part_2x2( T, ...
		0, 0, 'FLA_TL' );

    
	while ( size( ATL, 2 ) < size( A, 2  ) && size(ATL, 1) < size( A, 1 ) &&  num_iter > 0 )
		% Perform any column pivoting before partitioning
        if pivot_true == 1
		    pivot_idx = find_pivot(ABR);
            if pivot_idx ~= 1 % Pivot Index is relative to sL
			    base = size( ATL, 2 );
			    s([base + 1, pivot_idx + base]) = s([pivot_idx + base, base + 1]);
			    ABR(:, [1, pivot_idx]) = ABR(:, [pivot_idx, 1]);
			    ATR(:, [1, pivot_idx]) = ATR(:, [pivot_idx, 1]);
            end
        end

		[ A00,  a01,     A02,  ...
			a10t, alpha11, a12t, ...
			A20,  a21,     A22 ] = FLA_Repart_2x2_to_3x3( ATL, ATR, ...
														ABL, ABR, ...
														1, 1, 'FLA_BR' );
		[ T00,  ~,   T02,  ...
		  t10t, ~, t12t, ...
		  T20,  t21,   T22 ] = FLA_Repart_2x2_to_3x3( TTL, TTR, ...
														TBL, TBR, ...
														1, 1, 'FLA_BR' );
	
		%------------------------------------------------------------%
	
		[ alpha11, ...
			a21, tau11 ] = housev( alpha11, ...
									a21 );

        w12t = ( a12t + a21' * A22 )/ tau11;
		a12t = a12t - w12t;
		A22  = A22 - a21 * w12t;
		t01 = (a10t)' + A20'*a21;
	
		%------------------------------------------------------------%
	
		[ ATL, ATR, ...
			ABL, ABR ] = FLA_Cont_with_3x3_to_2x2( A00,  a01,     A02,  ...
												a10t, alpha11, a12t, ...
												A20,  a21,     A22, ...
												'FLA_TL' );
	


		[ TTL, TTR, ...
		  TBL, TBR ] = FLA_Cont_with_3x3_to_2x2( T00, t01, T02,  ...
												 t10t, tau11, t12t, ...
												 T20,  t21, T22, ...
												'FLA_TL' );
        num_iter = num_iter - 1;

	end
  
	A_out = [ ATL, ATR;
			  ABL, ABR ];
		  
	% t_out = [ tT,; 
	% 		  tB ];
  
	T_out = [TTL, TTR,; 
			 TBL, TBR];


  return
  
function [ idx ] = find_pivot(A)
	[~, idx ] = max(vecnorm(A, 2, 1));
return