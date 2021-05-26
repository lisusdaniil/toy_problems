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
Q = vel_var;
R = range_var;

% Define initial values
x_i = sqrt(max(r_dat.^2 - ell^2, 0.0));      % Take states to just be what you measured

% Get num of states
K = length(r_dat);

% Construct W matrix
Q_copy = repmat({Q\eye(size(Q))}, 1, K-1);
R_copy = repmat({R\eye(size(R))}, 1, K);
W = blkdiag(P_check_0\eye(size(P_check_0)), Q_copy{:}, R_copy{:});

% Start iterations
H = make_H(x_i, K, length(Q), ell);
e = make_e(x_i, K, ell, x_check_0, v_dat, r_dat, dt);
cost = 0.5*(e.')*e;
A = H'*W*H;
g = H'*W*e;
iter_num = 0;
nu = 2;
mu = 1e-6*max(diag(A));
found = false;
while ~found && iter_num < max_iter
    % Compute search direction
    del_x = -(A + mu*eye(size(A)))\g;
    
    % Check if step change too small
    if norm(del_x) < tol_step
        found = true;
    else
        % Update estimate
        x_new = x_i + del_x;

        % Get new error
        e_new = make_e(x_new, K, ell, x_check_0, v_dat, r_dat, dt);
        cost_new = 0.5*(e_new.')*e_new;

        gain_ratio = (cost - cost_new)/(0.5*del_x.'*(mu*del_x - g));
        if gain_ratio > 0       
            % Step accepted
            x_i = x_new;
            mu = mu*max(1/3,1 - (2*gain_ratio -1)^3);
            nu = 2;
            H = make_H(x_i, K, length(Q), ell);
            e = make_e(x_i, K, ell, x_check_0, v_dat, r_dat, dt);
            A = H'*W*H;
            g = H'*W*e;
            cost = 0.5*(e.')*e;
        else
            % Step rejected
            mu = mu*nu;
            nu = 2*nu;
        end
    end
    iter_num = iter_num + 1;
    
    if norm(g, inf) < tol_grad
       found = true; 
    end
end

% Evaluate x_hat
[m, n] = size(H);
Sigma = A\eye(n);
x_batch = x_i - Sigma*g;

%% Jacobian and Error construction
function H = make_H(x_i, K, len_Q, ell)
    H_top = eye(K*len_Q);
    for ii = 1:K-1
        H_top(ii*len_Q+1:(ii+1)*len_Q, 1+(ii-1)*len_Q:(ii)*len_Q) = -del_f_del_x_kmin(x_i(ii));
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
