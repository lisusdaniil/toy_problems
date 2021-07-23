% AUTHOR:       Daniil Lisus (260669654)
% TOPIC:        Nonlinear Batch Estimation
% MODULE:       Nonlinear Batch Estimation
% DESCRIPTION:  

% Choose which data to use
if plot_ode
    r_dat = r_ode_cor;
    v_dat = v_ode_cor;
else
    r_dat = r_sin_cor;
    v_dat = v_sin_cor;
end

% Define static matricies
A = 1;
Q = vel_var;
R = range_var;

% Get num of states to consider at a time
K = win_size;
K_m = win_size-marg_size;
%K = length(r_dat);

% Initialization
x_batch = zeros(length(r_dat),1);   % Initialize final estimate for plot
Sigma = zeros(length(r_dat));       % Initialize final cov estimate for plot

Q_copy = repmat({Q\eye(size(Q))}, 1, K-1);      % This doesnt change
R_copy = repmat({R\eye(size(R))}, 1, K);        % This doesnt change
Q_copy_m = repmat({Q\eye(size(Q))}, 1, ((K-K_m+1))-1);      % This doesnt change
R_copy_m = repmat({R\eye(size(R))}, 1, (K-K_m+1));        % This doesnt change

win_pop = 1;    % We initialize window with our initial guess
curr_win = zeros(K, 1);
curr_win(1) = x_check_0;

% Loop through each timestep
for t_idx = 2:length(T)
    % Check if window population is full
    if win_pop ~= (K-1)
        % initialize new window member from the range measurement
        win_pop = win_pop + 1;
        curr_win(win_pop) = sqrt(max(r_dat(t_idx,2).^2 - ell^2, 0.0));
        
    else
        % Construct W matrix using most up-to-date P_check_0
        W = blkdiag(P_check_0\eye(size(P_check_0)), Q_copy{:}, R_copy{:});
        
        % Populate last window member
        win_pop = win_pop + 1;
        curr_win(win_pop) = sqrt(max(r_dat(t_idx,2).^2 - ell^2, 0.0));
        
        % We get here if window is now populated
        % First, perform batch estimate on current window
        
        % Start iterations
        x_i = curr_win;
        iter_num = 0;
        del_x = 10*tol_step;
 
        % Time period here is from t_idx-(K-1) to t_idx
        win_t = (t_idx-(K-1)):t_idx;
        
        while norm(del_x) > tol_step && iter_num < max_iter
            % Construct H matrix
            % Construct H matrix
            H = make_H(x_i, K, length(A), ell);

            % Get error
            e = make_e(x_i, K, ell, x_check_0, v_dat(win_t,2), r_dat(win_t,2), dt);

            % Compute Gauss-Newton
            del_x = -(H'*W*H)\H'*W*e;

            % Update estimate
            alpha = 0.01;                % TODO: Implement line search algorithm
            x_i = x_i + alpha * del_x;
            iter_num = iter_num + 1;
        end

        % Evaluate x_hat
        Sigma(win_t,win_t) = (H'*W*H)\eye(size(H,2));
        x_batch(win_t) = x_i - Sigma(win_t,win_t)*H'*W*e;
       	
        % Second, marginalize out the respective states
        % Construct W_m matrix using most up-to-date P_check_0
        W_m = blkdiag(P_check_0\eye(size(P_check_0)), Q_copy_m{:}, R_copy_m{:});
        
        % Construct e_m from start of window to final marginalized state
        win_t_m = win_t(1:(K-K_m+1));
        e_m = make_e(x_batch(win_t_m), (K-K_m+1), ell, x_check_0, v_dat(win_t_m,2), r_dat(win_t_m,2), dt);
        
        % Construct H_m
        H_m = make_H(x_batch(win_t_m), (K-K_m+1), length(A), ell);
        
        % Isolate mu_m and Sigma_m
        Sigma_0m = (H_m'*W_m*H_m)\eye(size(H_m,2));
        Sigma_m = Sigma_0m(end,end);
        
        mu_0m = x_batch(win_t_m) - Sigma_0m*H_m'*W_m*e_m;
        mu_m = mu_0m(end);
        
        % Update new initial guess
        P_check_0 = Sigma_m;
        x_check_0 = mu_m;
        
        % Shift window
        curr_win = zeros(K, 1);
        curr_win(1) = x_check_0;
        curr_win(2:K_m) = x_i((K- K_m + 2):K);
        win_pop = K_m;
    end
end

%% Jacobian and Error construction
function H = make_H(x_i, K, len_A, ell)
    H_top = eye(K*len_A);
    for ii = 1:K-1
        H_top(ii*len_A+1:(ii+1)*len_A, 1+(ii-1)*len_A:(ii)*len_A) = -del_f_del_x_kmin(x_i(ii));
    end
    H_bot = -diag(del_g_del_x(x_i, ell));
    H = [H_top; H_bot];
end

function e = make_e(x_i, K, ell, x_check_0, v_dat, r_dat, dt)
    e = zeros(2*K,1);
    e(1) = x_i(1) - x_check_0;
    e(2:K) = x_i(2:K) - f_kmin(x_i(1:end-1), v_dat(1:end-1), dt);
    e(K+1:end) = r_dat - g_k(x_i, ell);
end

%% Process/Measurement models
function res = f_kmin(x_i, v_dat, dt)
    res = x_i + dt*v_dat;
end

function res = g_k(x_i, ell)
    res = sqrt(x_i.^2 + ell^2);
end

%% Partial functions
function res = del_f_del_x_kmin(x_kmin)
    res = 1;
end

% Compute full partial since we dont have to do weird off-diagonal stuff
function res = del_g_del_x(x_i, ell)
    res = x_i./sqrt(x_i.^2 + ell^2);
end
