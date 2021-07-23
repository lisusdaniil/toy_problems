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

% Initialize Q and R
Q_copy = repmat({Q\eye(size(Q))}, 1, K-1);      % This doesnt change
R_copy = repmat({R\eye(size(R))}, 1, K);        % This doesnt change
Q_copy_m = repmat({Q\eye(size(Q))}, 1, ((K-K_m+1))-1);      % This doesnt change
R_copy_m = repmat({R\eye(size(R))}, 1, (K-K_m+1));        % This doesnt change

% Initialization
x_batch = zeros(length(Q), length(test_path.r_meas));   % Initialize final estimate for plot
x_true = zeros(length(Q), length(test_path.r_meas));
err_batch = zeros(length(Q), length(test_path.r_meas));
sig_batch = zeros(length(Q), length(test_path.r_meas));

% Set up sliding window stuff
win_pop = 1;    % We initialize window with our initial guess
curr_win = zeros(4,4,K);
curr_win(:,:,1) = X_check_0;
%curr_win(:,:,1) = test_path.X_SE3(:,:,1)

% Loop through each timestep
blah = 0;
for t_idx = 2:length(T)
    
    % Check if window population is full
    if win_pop ~= (K-1)
        win_pop = win_pop + 1;
        % Dead reckon from previous state
        curr_win(:,:,win_pop) = curr_win(:,:,win_pop-1)*test_path.Xi_SE3(:,:,t_idx);
        %curr_win(:,:,win_pop) = test_path.X_SE3(:,:,t_idx);
    else
        % Construct W matrix using most up-to-date P_check_0
        W = blkdiag(P_check_0\eye(size(P_check_0)), Q_copy{:}, R_copy{:});
        
        % Populate last window member
        win_pop = win_pop + 1;
        curr_win(:,:,win_pop) = curr_win(:,:,win_pop-1)*test_path.Xi_SE3(:,:,t_idx);
        %curr_win(:,:,win_pop) = test_path.X_SE3(:,:,t_idx);
        
        % We get here if window is now populated
        % First, perform batch estimate on current window

        % Start iterations
        X_i = curr_win;
        iter_num = 0;
        
        % Time period for the window
        win_t = (t_idx-(K-1)):t_idx;
        
        % Extract relevant measurement matrices
        F_2 = [];
        F_1 = [];
        G_1 = [];
        for ii = 1:K
            F_2(:,:,ii) = eye(size(Q));
            F_1(:,:,ii) = -SE3.adjoint((test_path.Xi_SE3(:,:,win_t(ii)))\eye(4));
            G_1(:,:,ii) = [zeros(3), eye(3)];
        end
        
        % Set up LM
        [H, e] = make_H_e(X_i, K, length(Q), X_check_0, F_1, ...
            F_2, G_1, test_path.Xi_SE3(:,:,win_t),...
            test_path.r_meas(:,win_t));
        cost = 0.5*(e.')*e;
        A = H'*W*H;
        g = H'*W*e;
        nu = 1;
        mu = 0;
        found = false;
        
        while (~found && iter_num < max_iter)
            if iter_num == LM_iter
               mu = 1e+0*max(diag(A));
            end
            % Compute search direction
            del_xi = -(A + mu*eye(size(A)))\g;

            % Check if step change too small
            if norm(del_xi) < tol_step
                found = true;
            else
                % Update estimate
                X_new = zeros(size(X_i));
                for ii = 1:K
                    X_new(:,:,ii) = X_i(:,:,ii)*se3alg.expMap(-del_xi((ii-1)*length(Q)+1:ii*length(Q)));
                end

                % Get new error
                [H_new, e_new] = make_H_e(X_new, K, length(Q), X_check_0, F_1, F_2, G_1,...
                    test_path.Xi_SE3(:,:,win_t),...
                    test_path.r_meas(:,win_t));
                cost_new = 0.5*(e_new.')*e_new;

                gain_ratio = (cost - cost_new)/(0.5*del_xi.'*(mu*del_xi - g));
                if gain_ratio > 0 || iter_num < LM_iter  
                    % Step accepted
                    X_i = X_new;
                    mu = mu*max(1/3,1 - (2*gain_ratio -1)^3);
                    nu = 2;
                    H = H_new;
                    e = e_new;
                    cost = 0.5*(e.')*e;
                    A = H'*W*H;
                    g = H'*W*e;
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
        Sigma = (H'*W*H)\eye(size(H,2));
        for ii = 1:K
            x_batch(:,win_t(ii)) = se3alg.vee(SE3.logMap(X_i(:,:,ii)));
            x_true(:,win_t(ii)) = se3alg.vee(SE3.logMap(test_path.X_SE3(:,:,win_t(ii))));
            err_X = test_path.X_SE3(:,:,win_t(ii))\X_i(:,:,ii);
            err_batch(:,win_t(ii)) = se3alg.vee(SE3.logMap(err_X));
            sig_batch(:,win_t(ii)) = 3*diag(sqrt(Sigma((ii-1)*length(Q)+1:ii*length(Q),(ii-1)*length(Q)+1:ii*length(Q))));
        end

        % Second, marginalize out the respective states
        % Construct W_m matrix using most up-to-date P_check_0
        W_m = blkdiag(P_check_0\eye(size(P_check_0)), Q_copy_m{:}, R_copy_m{:});
        
        % Construct e_m and H_m from start of window to final marginalized state
        win_t_m = win_t(1:(K-K_m+1));
        
        [H_m, e_m] = make_H_e(X_i(:,:,1:(K-K_m+1)), (K-K_m+1), length(Q), X_check_0, F_1(:,:,1:(K-K_m+1)), ...
        F_2(:,:,1:(K-K_m+1)), G_1(:,:,1:(K-K_m+1)), test_path.Xi_SE3(:,:,win_t_m),...
        test_path.r_meas(:,win_t_m));
        
        % Isolate mu_m and Sigma_m
        iso_m_idx = (K-K_m)*length(Q)+1:(K-K_m+1)*length(Q);  % isolate desired indeces
        
        Sigma_0m = (H_m'*W_m*H_m)\eye(size(H_m,2));
        Sigma_m = Sigma_0m(iso_m_idx, iso_m_idx);
        
        % Get the perturbation on our state
        del_X_m = -Sigma_0m*H_m'*W_m*e_m;
        % Update X_m through perturbation
        mu_m = X_i(:,:,(K-K_m+1))*SE3.expMap(se3alg.wedge(del_X_m(iso_m_idx)));

        % Update new initial guess
        P_check_0 = Sigma_m;
        X_check_0 = mu_m;
        
        % Shift window
        curr_win = zeros(4,4,K);
        curr_win(:,:,1) = X_check_0;
        curr_win(:,:,2:K_m) = X_i(:,:,(K-K_m+2):K);
        win_pop = K_m;
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
