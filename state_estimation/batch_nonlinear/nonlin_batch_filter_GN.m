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
R = 10*range_var;

% Define initial values
x_i = sqrt(max(r_dat.^2 - ell^2, 0.0));      % Take states to just be what you measured

% Get num of states
K = length(r_dat);

% Construct W matrix
Q_copy = repmat({Q\eye(size(Q))}, 1, K-1);
R_copy = repmat({R\eye(size(R))}, 1, K);
W = blkdiag(P_check_0\eye(size(P_check_0)), Q_copy{:}, R_copy{:});

% Start iterations
iter_num = 0;
del_x = 10*tol_step;
while norm(del_x) > tol_step && iter_num < max_iter
    % Construct H matrix
    % Construct H matrix
    H = make_H(x_i, K, length(A), ell);
    
    % Get error
    e = make_e(x_i, K, ell, x_check_0, v_dat, r_dat, dt);
    
    % Compute Gauss-Newton
    del_x = -(H'*W*H)\H'*W*e;
    
    % Update estimate
    alpha = 0.01;                % TODO: Implement line search algorithm
    x_i = x_i + alpha * del_x;
    iter_num = iter_num + 1;
end

% Evaluate x_hat
[m, n] = size(H);
Sigma = (H'*W*H)\eye(n);
x_batch = x_i - Sigma*H'*W*e;

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
