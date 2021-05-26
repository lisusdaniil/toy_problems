% AUTHOR:       Daniil Lisus (260669654)
% TOPIC:        MLG Batch Estimation
% MODULE:       MLG Batch Estimator
% DESCRIPTION:  

% Get num of states to consider at a time
K = win_size;
K_m = win_size-marg_size;
%K = size(test_path.X_SE3,3);

% Define static matricies
Q = diag([gyro_x_std^2, gyro_y_std^2, gyro_z_std^2, vel_x_std^2, vel_y_std^2, vel_z_std^2]);
R = diag([GPS_x_std^2, GPS_y_std^2, GPS_z_std^2]);

Q = Q + eye(size(Q))*10^-8;
R = R + eye(size(R))*10^-8;


% Define the static H matrix elements
F_2 = [];
F_1 = [];
G_1 = [];
for ii = 1:K
    F_2(:,:,ii) = eye(size(Q));
    F_1(:,:,ii) = -SE3.adjoint((test_path.Xi_SE3(:,:,ii))\eye(4));
    G_1(:,:,ii) = [zeros(3), eye(3)];
end

% Initialization
x_batch = zeros(length(r_dat),1);   % Initialize final estimate for plot
Sigma = zeros(length(r_dat));       % Initialize final cov estimate for plot

Q_copy = repmat({Q\eye(size(Q))}, 1, K-1);      % This doesnt change
R_copy = repmat({R\eye(size(R))}, 1, K);        % This doesnt change
Q_copy_m = repmat({Q\eye(size(Q))}, 1, ((K-K_m+1))-1);      % This doesnt change
R_copy_m = repmat({R\eye(size(R))}, 1, (K-K_m+1));        % This doesnt change

win_pop = 0;
curr_win = zeros(4,4,K);

% Loop through each timestep
for t_idx = 1:length(T)
    t_k = T(t_idx);
    
    % Check if window population is full
    if win_pop ~= (K-1)
        % initialize new window member from the range measurement
        win_pop = win_pop + 1;
        curr_win(:,:, win_pop) = SE3.synthesize(zeros(3,1), test_path.r_meas(:,t_idx));
    else
        % Construct W matrix using most up-to-date P_check_0
        W = blkdiag(P_check_0\eye(size(P_check_0)), Q_copy{:}, R_copy{:});
        
        % Populate last window member
        win_pop = win_pop + 1;
        curr_win(:,:,win_pop) = SE3.synthesize(zeros(3,1), test_path.r_meas(:,t_idx));
        
        % We get here if window is now populated
        % First, perform batch estimate on current window

        % Start iterations
        X_i = curr_win;
        iter_num = 0;
        del_xi = 10*tol_step;
        
        % Time period here is from t_idx-(K-1) to t_idx
        win_t = (t_idx-(K-1)):t_idx;
        
        while norm(del_xi) > tol_step && iter_num < max_iter
            % Get Jacobian and error
            [H, e] = make_H_e(X_i, K, length(Q), X_check_0, F_1(:,:,win_t), ...
            F_2(:,:,win_t), G_1(:,:,win_t), test_path.Xi_SE3(:,:,win_t),...
            test_path.r_meas(:,win_t));

            % Compute Gauss-Newton
            del_xi = -(H'*W*H)\H'*W*e;

            % Update estimate
            alpha = 0.1;                % TODO: Implement line search algorithm
            for ii = win_t
                X_i(:,:,ii) = X_i(:,:,ii)*se3alg.expMap(-alpha * del_xi((ii-1)*length(Q)+1:ii*length(Q)));
            end
            iter_num = iter_num + 1;
        end

        % Evaluate x_hat
        [m, n] = size(H);
        Sigma = (H'*W*H)\eye(n);
        x_batch = [];
        x_true = [];
        err_X = [];
        err_batch = [];
        sig_batch = [];
        for ii = 1:K
            x_batch(:,ii) = se3alg.vee(SE3.logMap(X_i(:,:,ii)));
            x_true(:,ii) = se3alg.vee(SE3.logMap(test_path.X_SE3(:,:,ii)));
            err_X = test_path.X_SE3(:,:,ii)\X_i(:,:,ii);
            err_batch(:,ii) = se3alg.vee(SE3.logMap(err_X));
            sig_batch(:,ii) = 3*diag(sqrt(Sigma((ii-1)*length(Q)+1:ii*length(Q),(ii-1)*length(Q)+1:ii*length(Q))));
        end

    end
end
        
%% Function to evaluate Jacobians and errors
function [H, e] = make_H_e(X_i, K, len_Q, X_check_0, F_1, F_2, G_1, Xi_SE3, r_meas)
    % Construct H matrix
    H_top = zeros((K)*len_Q);
    H_bot = zeros((K)*(size(X_i,1)-1), (K)*len_Q);
    
    % Construct errors
    e_0 = se3alg.vee(SE3.logMap(X_i(:,:,1)\X_check_0));
    e_u = zeros([(K-1)*len_Q,1]);
    e_y = zeros([(K)*(size(X_i,1)-1),1]);
    
    for ii = 1:K
        % Compose top of H
        H_top((ii-1)*(len_Q)+1:ii*(len_Q), (ii-1)*(len_Q)+1:ii*(len_Q)) = F_2(:,:,ii);
        if ii < K
            H_top(ii*len_Q+1:(ii+1)*len_Q, 1+(ii-1)*len_Q:(ii)*len_Q) = F_1(:,:,ii);
        end

        % Compose bot of H
        % Get most up to date G_1
        H_bot((ii-1)*(size(X_i,1)-1)+1:ii*(size(X_i,1)-1), (ii-1)*(len_Q)+1:ii*(len_Q)) = G_1(:,:,ii);
        
        % Get errors
        if ii > 1
            E_u = X_i(:,:,ii)\(X_i(:,:,ii-1)*Xi_SE3(:,:,ii-1));
            e_u((ii-2)*(len_Q)+1:(ii-1)*(len_Q)) = se3alg.vee(SE3.logMap(E_u));
        end
        E_y = X_i(:,:,ii)\([r_meas(:,ii); 1.0] - X_i(:,:,ii)*[zeros(3,1); 1]);
        e_y((ii-1)*(size(X_i,1)-1)+1:ii*(size(X_i,1)-1)) = E_y(1:3);
    end
    
    H = [H_top; H_bot];
    e = [e_0; e_u; e_y];
end
