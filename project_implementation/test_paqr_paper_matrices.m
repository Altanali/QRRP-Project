addpath('paqr matlab source/');
matrixNameIdx = 1;
matrixFunIdx  = 2;
matrixInfo = { ...
    {'rand'         @(n) rand(n,n) }
    {'vander'       @(n) vander(linspace(0, 1, n)) }
    {'baart'        @baart       }
    {'break1'       @break1      }
    {'break9'       @break9      }
    {'deriv2'       @deriv2      }
    {'devil'        @devil       }
    {'exponential'  @exponential }
    {'foxgood'      @foxgood     }
    {'gravity'      @gravity     }
    {'heat'         @heat        }
    {'phillips'     @phillips    }
    {'random'       @random      }
    {'shaw'         @shaw        }
    {'spikes'       @spikes      }
    {'stewart'      @stewart     }
    {'ursell'       @ursell      } 
    {'wing'         @wing        }
};

n = 1000;
num_iter = 4;
headers = {"matrix name", "paqr", "hqrrp", "pa_hqrrp", "hqr", "hqrp"};

% Initialize result matrices
time_results = cell(length(matrixInfo)+1, length(headers));
error_results = cell(length(matrixInfo)+1, length(headers));

time_results(2:end, 2:end) = {0};
error_results(2:end, 2:end) = {0};

% Add headers to the result matrices
time_results(1, :) = headers;
error_results(1, :) = headers;

% Loop through each matrix and run tests
for i = 1:length(matrixInfo)
    matrixName = matrixInfo{i}{matrixNameIdx};
    matrixFun  = matrixInfo{i}{matrixFunIdx};

    
    % Store the matrix name
    time_results{i+1, 1} = matrixName;
    error_results{i+1, 1} = matrixName;
    fprintf("%s\n", matrixName);
    
    for iter = 1:num_iter
        A = matrixFun(n);
        x = rand(n, 1);
        
        % Test PAQR
        [error, time] = test_qr_algo(A, x, "paqr");
        time_results{i+1, 2} = time_results{i+1, 2} + time;
        error_results{i+1, 2} = error_results{i+1, 2}+ error;
    
        % Test HQRRP
        [error, time] = test_qr_algo(A, x, "hqrrp");
        time_results{i+1, 3} = time_results{i+1, 3} + time;
        error_results{i+1, 3} = error_results{i+1, 3} + error;
    
        % Test PA_HQRRP
        [error, time] = test_qr_algo(A, x, "pa_hqrrp");
        time_results{i+1, 4} = time_results{i+1, 4} + time;
        error_results{i+1, 4} = error_results{i+1, 4} + error;
    
        % Test HQR
        [error, time] = test_qr_algo(A, x, "hqr");
        time_results{i+1, 5} = time_results{i+1, 5} + time;
        error_results{i+1, 5} = error_results{i+1, 5} + error;
    
        % Test HQRP
        [error, time] = test_qr_algo(A, x, "hqrp");
        time_results{i+1, 6} = time_results{i+1, 6} + time;
        error_results{i+1, 6} = error_results{i+1, 6} + error;
    end
    for algo_idx = 2:6
        time_results{i+1, algo_idx} = time_results{i+1, algo_idx} ./num_iter;
        error_results{i+1, algo_idx} = error_results{i+1, algo_idx} ./num_iter;
    end
end

% Write the results to CSV files
writecell(time_results, 'project_implementation/data/paqr_matrices_time_results.csv');
writecell(error_results, 'project_implementation/data/paqr_matrices_error_results.csv');