% AUTHOR:       Daniil Lisus (260669654)
% TOPIC:        MLG Batch Estimation
% MODULE:       MLG Batch Estimator
% DESCRIPTION:  

% Get num of states
K = size(test_path.X_SE3,3);

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

% Define initial values
X_i = [];
for ii = 1:K
    X_i(:,:,ii) = SE3.synthesize(zeros(3,1), test_path.r_meas(:,ii));
    %X_i(:,:,ii) = test_path.X_SE3(:,:,ii);
end


% Construct W matrix
Q_copy = repmat({Q\eye(size(Q))}, 1, K-1);
R_copy = repmat({R\eye(size(R))}, 1, K);
W = blkdiag(P_check_0\eye(size(P_check_0)), Q_copy{:}, R_copy{:});

% Start iterations
[H, e] = make_H_e(X_i, K, length(Q), X_check_0, F_1, F_2, G_1, test_path);
cost = 0.5*(e.')*e;
A = H'*W*H;
g = H'*W*e;
iter_num = 0;
nu = 2;
mu = 0;
found = false;
while ~found && iter_num < max_iter
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
        [H_new, e_new] = make_H_e(X_new, K, length(Q), X_check_0, F_1, F_2, G_1, test_path);
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

%% Function to evaluate Jacobians and errors
function [H, e] = make_H_e(X_i, K, len_Q, X_check_0, F_1, F_2, G_1, test_path)
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
            E_u = X_i(:,:,ii)\(X_i(:,:,ii-1)*test_path.Xi_SE3(:,:,ii-1));
            e_u((ii-2)*(len_Q)+1:(ii-1)*(len_Q)) = se3alg.vee(SE3.logMap(E_u));
        end
        E_y = X_i(:,:,ii)\([test_path.r_meas(:,ii); 1.0] - X_i(:,:,ii)*[zeros(3,1); 1]);
        e_y((ii-1)*(size(X_i,1)-1)+1:ii*(size(X_i,1)-1)) = E_y(1:3);
    end
    
    H = [H_top; H_bot];
    e = [e_0; e_u; e_y];
end
