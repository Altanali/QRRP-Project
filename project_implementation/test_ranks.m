results = struct('test_paqr', [], 'test_hqrrp', [], 'test_pa_hqrrp', [], 'test_hqr', [], 'test_hqrp', []);
m = 1000;
num_iter = 2;
for rk = 100:100:m
    fprintf("=========== Testing Rank: %d ===========\n", rk);

    paqr_time = 0;
    paqr_error = 0;

    pa_hqrrp_time = 0;
    pa_hqrrp_error = 0;

    hqrrp_time = 0;
    hqrrp_error = 0;

    hqrp_time = 0;
    hqrp_error = 0;

    hqr_time = 0;
    hqr_error = 0;

    for iter=1:num_iter
        fprintf("------ (rk=%d) Iteration %d ------\n", rk, iter);
        % Generate random matrices and vectors
        A = zeros(m, m); 
        %Create a rank rk matrix via a sum of rk outer products. 
        for i = 1:rk
            x_temp = rand(m, 1);
            A = A + x_temp * x_temp';
        end
        x = rand(m, 1);
        % Run each test function and store the results
        [error, time] = test_qr_algo(A, x, "paqr");
        paqr_time = paqr_time + time;
        paqr_error = paqr_error +  error;

        [error, time] = test_qr_algo(A, x, "hqrrp");    
        hqrrp_time = hqrrp_time + time;
        hqrrp_error = hqrrp_error +  error;

        [error, time] = test_qr_algo(A, x, "pa_hqrrp");    
        pa_hqrrp_time = pa_hqrrp_time + time;
        pa_hqrrp_error = pa_hqrrp_error +  error;

        [error, time]  = test_qr_algo(A, x, "hqr");   
        hqr_time = hqr_time + time;
        hqr_error = hqr_error +  error;

        [error, time]  = test_qr_algo(A, x, "hqrp");
        hqrp_time = hqrp_time + time;
        hqrp_error = hqrp_error +  error;
    end
    % Store the errors and times in arrays
    errors = [paqr_error, hqrrp_error, pa_hqrrp_error, hqr_error, hqrp_error];
    times = [paqr_time, hqrrp_time, pa_hqrrp_time, hqr_time, hqrp_time];
    
    % Divide all errors and times by num_iter
    errors = errors ./ num_iter;
    times = times ./ num_iter;
    
    % Append the results to the corresponding fields in the results struct
    results.test_paqr = [results.test_paqr; rk, times(1), errors(1)];
    results.test_hqrrp = [results.test_hqrrp; rk, times(2), errors(2)];
    results.test_pa_hqrrp = [results.test_pa_hqrrp; rk, times(3), errors(3)];
    results.test_hqr = [results.test_hqr; rk, times(4), errors(4)];
    results.test_hqrp = [results.test_hqrp; rk, times(5), errors(5)];
end

% Save the results to a .mat file for later use
save('project_implementation/results.mat', 'results');
%%