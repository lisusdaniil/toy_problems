% AUTHOR:       Daniil Lisus (260669654)
% TOPIC:        Nonlinear IMU Preintegration Estimation
% MODULE:       Nonlinear IMU Preintegration GN Estimation
% DESCRIPTION:  

% Choose which data to use
if plot_ode
    r_dat = r_ode_cor;
    v_dat = v_ode_cor;
else
    r_dat = r_sin_cor;
    v_dat = v_sin_cor;
end

% Downsample extero
r_dat = r_dat(1:downsample_freq:end);

% Define static matricies
Q = vel_var;
R = range_var;

% Get num of states
K_int = length(v_dat);
K_ext = length(r_dat);
if IMU_preint
    [v_input, Q_copy] = preint_vel(v_dat, Q, dt, K_int, K_ext);
    K_int = length(v_input);    % Rewrite new number of intero inputs
    T_batch = T(1:downsample_freq:end);
else
    % Define standard odom input and uncertainty
    v_input = dt*v_dat;
    Q_copy = repmat({Q\eye(size(Q))}, 1, K_int-1);
    T_batch = T;
end

% Define initial values
x_i = ones(K_int,1);      % Take states to just be what you measured

% Construct W matrix
R_copy = repmat({R\eye(size(R))}, 1, K_ext);
W = blkdiag(P_check_0\eye(size(P_check_0)), Q_copy{:}, R_copy{:});

% Start iterations
iter_num = 0;
del_x = 10*tol_step;
while norm(del_x) > tol_step && iter_num < max_iter
    % Construct H matrix
    % Construct H matrix
    H = make_H(x_i, K_int, K_ext, ell);
    
    % Get error
    e = make_e(x_i, K_int, K_ext, ell, x_check_0, v_input, r_dat);
    
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
function H = make_H(x_i, K_int, K_ext, ell)
    H_top = eye(K_int);
    for ii = 1:K_int-1
        H_top(ii+1:(ii+1), 1+(ii-1):(ii)) = -del_f_del_x_kmin(x_i(ii));
    end
    H_bot = zeros(K_ext, K_int);
    for ii = 1:K_ext
        range_idx = 1+(ii-1)*(K_int-1)/(K_ext-1);
        H_bot(ii, range_idx) = -del_g_del_x(x_i(range_idx), ell);
    end
    H = [H_top; H_bot];
end

function e = make_e(x_i, K_int, K_ext, ell, x_check_0, v_input, r_dat)
    e = zeros(K_int+K_ext,1);
    e(1) = x_i(1) - x_check_0;
    e(2:K_int) = x_i(2:K_int) - f_kmin(x_i(1:end-1), v_input(1:end-1));
    
    for ii = 1:K_ext
        range_idx = 1+(ii-1)*(K_int-1)/(K_ext-1);
        e(K_int+ii) = r_dat(ii) - g_k(x_i(range_idx), ell);
    end
end

%% Process/Measurement models
function res = f_kmin(x_i, v_input)
    res = x_i + v_input;
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

%% Vel Preintegration
function [v_input, Sigmas] = preint_vel(v_dat, Q, dt, K_int, K_ext)
    % Initialize final matrices
    v_input = zeros(K_ext,1);
    Sigmas = cell(1, K_ext);
    
    for ii = 1:K_ext-1
        % Find indeces from which to which to preintegrate
        idx_k = 1+(ii-1)*(K_int-1)/(K_ext-1);
        idx_k1 = (ii)*(K_int-1)/(K_ext-1);
        
        v_tmp = 0;
        Sig_tmp = 0;
        F = 0;
        for jj = idx_k:idx_k1
            v_tmp = v_tmp + v_dat(jj)*dt;
            Sig_tmp = Sig_tmp + (1 + dt*F)^2*Sig_tmp + dt*Q;
        end
        
        v_input(ii) = v_tmp;
        Sigmas{ii} = Sig_tmp;
    end
end
