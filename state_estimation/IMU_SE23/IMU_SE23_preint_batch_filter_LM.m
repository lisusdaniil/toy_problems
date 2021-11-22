% AUTHOR:       Daniil Lisus (260669654)
% TOPIC:        SE2(3) with bias IMU Preintegration Simulated Data
% MODULE:       LM SWF Estimation
% DESCRIPTION:  

%% Unit tests
%proc_err_result = process_error_unit_test();

%%
% Isolate measurements
gyr_dat = sensorData.meas_gyro;
a_dat = sensorData.meas_accel;
v_lin_dat = sensorData.meas_odom_lin;
v_ang_dat = sensorData.meas_odom_ang;
r_dat = sensorData.meas_position;

%gyr_dat = states_imu.omega_bt_b;
%a_dat = states_imu.f_b_tan;
%v_lin_dat = states_odom.v_btt_b;
%v_ang_dat = states_odom.omega_bt_b;
%r_dat = gt_states.r_bt_t;

% Initialize from settings
X_check_0 = X_check_0_init;
P_check_0 = P_check_0_init;

% Get num of states to consider at a time
K = win_size;
K_m = win_size-marg_size;

% Initialize uncertainty matricies
Q_r = diag([sigma_gyro^2, sigma_gyro^2, 5*sigma_gyro^2,...
            sigma_accel^2, sigma_accel^2, sigma_accel^2,...
            sigma_gyro_bias^2, sigma_gyro_bias^2, sigma_gyro_bias^2,...
            sigma_accel_bias^2, sigma_accel_bias^2, sigma_accel_bias^2]); % Continuous time PSD

% Define odometry uncertainty
odom_lin_ps_sigma = 0.01;    % [m / s], uncertainty on pseudo linear measurements
odom_ang_ps_sigma = 0.01;    % [rad / s], uncertainty on pseudo ang measurements
Q_o_r = 2*diag([odom_ang_ps_sigma^2, odom_ang_ps_sigma^2, sigma_odom_ang^2, sigma_odom_lin^2, odom_lin_ps_sigma^2, odom_lin_ps_sigma^2]);
% Define GPS measurement uncertainty
R = diag([sigma_position^2, sigma_position^2, 5*sigma_position^2]);

% Initialize W matrix
W_Q = P_check_0\eye(size(P_check_0));
W_Q_o = [];
W_R = R\eye(size(R));   % Initialize measurement W separately to add it on later
r_meas = r_dat(:,1);

% Initialization
num_xi = length(P_check_0);
if use_marg
    num_states = length(gt_states.time);    % Data gen takes care of IMU preint flag
else
    num_states = win_size;
end
dt_imu = sensorData.t_IMU(2);

x_batch = zeros(num_xi, num_states);
x_true = zeros(num_xi, num_states);
err_batch = zeros(num_xi, num_states);
sig_batch = zeros(num_xi, num_states);

% Set up sliding window stuff
win_pop = 1;    % We initialize window with our initial guess
curr_win = zeros([size(X_check_0),K]);
curr_win(:,:,1) = X_check_0;
del_X = [];
Del_t = [];
b_bar = [];
Sigma_del_X = [];
num_windows = 0;
del_X_o = [];
Sigma_del_X_o = [];

% Loop through each timestep index
for t_idx = 2:num_states
    % If IMU preint is on, then t_idx corresponds to extero time
    % Else, it corresponds to intero time
    win_pop = win_pop + 1;
    if IMU_preint
        % Preint all intero measurements between t_idx-1 to t_idx
        t_k = find_closest_t_idx(sensorData.t_IMU, sensorData.t_pos(t_idx-1));
        t_k_1 = find_closest_t_idx(sensorData.t_IMU, sensorData.t_pos(t_idx));
        t_int = t_k:(t_k_1-1);
        [del_X(:,:,win_pop-1), Sigma_del_X, Del_t(win_pop-1), bias_J(:,:,:,win_pop-1), b_bar(:,win_pop-1)] ...
            = preint_IMU(curr_win(:,:,win_pop-1), gyr_dat(:,t_int), a_dat(:,t_int), Q_r, dt_imu, error_deff);
        if use_odom
            [del_X_o(:,:,win_pop-1), Sigma_del_X_o] = preint_odom(v_lin_dat(:,t_int), v_ang_dat(:,t_int), Q_o_r, dt_imu, error_deff);
        end
    else
        del_C = so3alg.expMap(so3alg.wedge(gyr_dat(:,t_idx-1)*dt_imu));
        del_v = a_dat(:,idx-1)*dt_imu;
        del_r = 0.5*a_dat(:,idx-1)*dt_imu^2;
        del_b = zeros(3,1);
        % Assemble dead reckoned change in state
        del_X(:,:,t_idx-1) = SE23R3R3_synthesize(del_C, del_v, del_r, del_b, del_b);
        %Sigma_del_X = Q; Need to properly conver Q to Q_d
    end
    
    % Update W with new uncertainty from process model
    W_Q = blkdiag( W_Q, Sigma_del_X\eye(size(Sigma_del_X)) );
    % Update measurement part of W if new measurement available
    if use_odom
        W_Q_o = blkdiag( W_Q_o, Sigma_del_X_o\eye(size(Sigma_del_X_o)));
    end
    W_R = blkdiag( W_R, R\eye(size(R)) );
    r_meas = [r_meas, r_dat(:,t_idx)];
    
    % Check if window population is full
    if win_pop ~= K
        % Dead reckon from previous state
        curr_win(:,:,win_pop) = F(curr_win(:,:,win_pop-1), del_X(:,:,win_pop-1), Del_t(win_pop-1), g_a);
    else
        % Populate last window member
        curr_win(:,:,win_pop) = F(curr_win(:,:,win_pop-1), del_X(:,:,win_pop-1), Del_t(win_pop-1), g_a);
        %curr_win(:,:,win_pop) = test_path.X_SE3(:,:,t_idx);
        
        % We get here if window is now populated
        % First, perform batch estimate on current window

        % Start iterations
        X_i = curr_win;
        iter_num = 0;
        % Save bias_bar 
        del_X_og = del_X;
        
        % Time period for the window
        win_t = (t_idx-(K-1)):t_idx;
        
        % Form error Jacobian terms
        F_1 = [];
        F_2 = [];
        G_1 = [];
        G_o_1 = [];
        G_o_2 = [];
        [F_1, F_2, G_1] = make_err_Jacobians(X_i, del_X, Del_t, bias_J, b_bar, r_meas, error_deff);
        if use_odom
            [G_o_1, G_o_2] = make_odom_err_Jacobians(X_i, del_X_o, error_deff);
        end

        % Assemble full W
        W = blkdiag( W_Q, W_Q_o, W_R );
        
        % Set up LM
        H = make_H(F_1, F_2, G_1, G_o_1, G_o_2, use_odom);
        e = make_e(X_i, X_check_0, del_X, Del_t, bias_J, b_bar, g_a, r_meas, K, num_xi, del_X_o, use_odom, error_deff);
        cost = 0.5*(e.')*e;
        A = H'*W*H;
        g = H'*W*e;
        nu = 1;
        mu = 0;
        found = false;
        cost_idx = 1;

        while (~found && iter_num < max_iter)            
            if iter_num == LM_iter
               mu = 1e+0*max(diag(A));
            end
            % Compute search direction
            del_xi = -(A + mu*eye(size(A)))\g;
            
            % Display header
            if verbose
                header = ['Iter: ', num2str(iter_num), '   ||   ' ...
                          'Cost: ', num2str(cost), '   ||   '...
                          'Step: ', num2str(norm(del_xi))];

                disp(header)
            end

            % Check if step change too small
            if norm(del_xi) < tol_step
                found = true;
            else
                % Update estimate
                X_new = zeros(size(X_i));
                for ii = 1:K
                    del_xi_i = del_xi((ii-1)*num_xi+1:ii*num_xi);
                    if error_deff == "left"
                        X_new(:,:,ii) = X_i(:,:,ii) * se23xR3x2alg.expMap(se23xR3x2alg.wedge( -del_xi_i ));
                    else
                        X_new(:,:,ii) = se23xR3x2alg.expMap(se23xR3x2alg.wedge( -del_xi_i )) * X_i(:,:,ii);
                    end
                end
                
                F_1_new = [];
                F_2_new = [];
                G_1_new = [];
                G_o_1_new = [];
                G_o_2_new = [];
                [F_1_new, F_2_new, G_1_new] = make_err_Jacobians(X_new, del_X, Del_t, bias_J, b_bar, r_meas, error_deff);
                if use_odom
                    [G_o_1_new, G_o_2_new] = make_odom_err_Jacobians(X_new, del_X_o, error_deff);
                end
                
                H_new = make_H(F_1_new, F_2_new, G_1_new, G_o_1_new, G_o_2_new, use_odom);
                e_new = make_e(X_new, X_check_0, del_X, Del_t, bias_J, b_bar, g_a, r_meas, K, num_xi, del_X_o, use_odom, error_deff);
                
                cost_new = 0.5*(e_new.')*e_new;

                gain_ratio = (cost - cost_new)/(0.5*del_xi.'*(mu*del_xi - g));
                if iter_num == LM_iter
                   display("LM started") 
                end
                if gain_ratio > 0 || iter_num < LM_iter  
                    % Step accepted
                    X_i = X_new;
                    
                    % Propagate bias
                    %del_X = prop_bias(del_X_og, bias_J, del_b);
                    % Form error Jacobian terms
                    [F_1, F_2, G_1] = make_err_Jacobians(X_i, del_X, Del_t, bias_J, b_bar, r_meas, error_deff);
                    if use_odom
                        [G_o_1, G_o_2] = make_odom_err_Jacobians(X_i, del_X_o, error_deff);
                    end
                    
                    mu = mu*max(1/3,1 - (2*gain_ratio -1)^3);
                    nu = 2;
                    H = H_new;
                    e = e_new;
                    cost = 0.5*(e.')*e;
                    cost_save(cost_idx) = cost;
                    A = H'*W*H;
                    g = H'*W*e;
                    cost_idx = cost_idx + 1;
                else
                    % Step rejected
                    mu = mu*nu;
                    nu = 2*nu;
                    disp("Rejected Step")
                end
            end
            iter_num = iter_num + 1;

            if norm(g, inf) < tol_grad
               found = true; 
            end
        end
        
        %% Evaluate x_hat
        Sigma = (H'*W*H)\eye(size(H,2));
        for ii = 1:K
            x_batch(:,win_t(ii)) = se23xR3x2alg.vee( SE23xR3x2.logMap( X_i(:,:,ii) ) );
            x_true(:,win_t(ii)) = se23xR3x2alg.vee( SE23xR3x2.logMap( gt_SE23R3R3(:,:,win_t(ii)) ) );
            err_X = gt_SE23R3R3(:,:,win_t(ii)) \ X_i(:,:,ii);
            err1 = se23xR3x2alg.vee(SE23xR3x2.logMap(gt_SE23R3R3(:,:,win_t(ii)) \ X_i(:,:,ii)));
            err2 = se23xR3x2alg.vee(SE23xR3x2.logMap( X_i(:,:,ii) / gt_SE23R3R3(:,:,win_t(ii))));
            
            err_batch(:,win_t(ii)) = se23xR3x2alg.vee(SE23xR3x2.logMap(err_X));
            Sigma_ii = Sigma((ii-1)*num_xi+1:ii*num_xi,(ii-1)*num_xi+1:ii*num_xi);
            if error_deff == "right"
                % If right invariant error, transform uncertainty to left
                % for plotting
                Sigma_ii = inv(SE23xR3x2.adjoint(X_i(:,:,ii))) * Sigma_ii * inv(SE23xR3x2.adjoint(X_i(:,:,ii))');
            end
            sig_batch(:,win_t(ii)) = 3*diag(sqrt(Sigma_ii));
            Sig_batch(:,:,win_t(ii)) = Sigma_ii;
        end
        
        if ~use_marg
           break; 
        end
        
        % Second, marginalize out the respective states
        % Construct W_m matrix using most up-to-date P_check_0
        W_Q_marg_span = 1:(K-K_m+1)*num_xi;
        W_R_marg_span = 1:(K-K_m+1)*length(R);
        W_Q_o_marg_span = 1:(K-K_m)*length(Sigma_del_X_o);
        if use_odom
            W_m = blkdiag(W_Q(W_Q_marg_span,W_Q_marg_span), W_Q_o(W_Q_o_marg_span,W_Q_o_marg_span), W_R(W_R_marg_span,W_R_marg_span));
        else
            W_m = blkdiag(W_Q(W_Q_marg_span,W_Q_marg_span), W_R(W_R_marg_span,W_R_marg_span));
        end
        
        % Construct e_m and H_m from start of window to final marginalized state
        X_i_m = X_i(:,:,1:(K-K_m+1));
        del_X_m = del_X(:,:,1:(K-K_m));
        Del_t_m = Del_t(1:(K-K_m));
        bias_J_m = bias_J(:,:,:,1:(K-K_m));
        b_bar_m = b_bar(:,1:(K-K_m));
        F_1_m = F_1(:,:,1:(K-K_m));     % One less than in F_2_m since F_0
        F_2_m = F_2(:,:,1:(K-K_m+1));
        G_1_m = G_1(:,:,1:(K-K_m+1));
        r_meas_m = r_meas(:,1:(K-K_m+1));
        if use_odom
            G_o_1_m = G_o_1(:,:,1:(K-K_m));
            G_o_2_m = G_o_2(:,:,1:(K-K_m));
            del_X_o_m = del_X_o(:,:,1:(K-K_m));
        else
            G_o_1_m = [];
            G_o_2_m = [];
            del_X_o_m = [];
        end
        
        H_m = make_H(F_1_m, F_2_m, G_1_m, G_o_1_m, G_o_2_m, use_odom);
        e_m = make_e(X_i_m, X_check_0, del_X_m, Del_t_m, bias_J_m, b_bar_m, g_a, r_meas_m, (K-K_m+1), num_xi, del_X_o_m, use_odom, error_deff);
        
        % Isolate mu_m and Sigma_m
        iso_m_idx = (K-K_m)*num_xi+1:(K-K_m+1)*num_xi;  % isolate desired indices
        
        Sigma_0m = (H_m'*W_m*H_m)\eye(size(H_m,2));
        Sigma_m = Sigma_0m(iso_m_idx, iso_m_idx);
        
        % Get the perturbation on our state
        del_xi_m = -Sigma_0m*H_m'*W_m*e_m;
        % Update X_m through perturbation
        if error_deff == "left"
            mu_m = X_i_m(:,:,end) * se23xR3x2alg.expMap(se23xR3x2alg.wedge(-del_xi_m(iso_m_idx)));
        else
            mu_m = se23xR3x2alg.expMap(se23xR3x2alg.wedge(-del_xi_m(iso_m_idx))) *  X_i_m(:,:,end);
        end

        % Update new initial guess
        P_check_0 = Sigma_m;
        X_check_0 = mu_m;
        
        % Shift window
        curr_win = zeros(length(X_check_0),length(X_check_0),K);
        curr_win(:,:,1) = X_check_0;
        curr_win(:,:,2:K_m) = X_i(:,:,(K-K_m+2):K);
        del_X =  del_X(:,:,(K-K_m+1):K-1);
        Del_t = Del_t((K-K_m+1):K-1);
        bias_J = bias_J(:,:,:, (K-K_m+1):K-1);
        b_bar = b_bar(:, (K-K_m+1):K-1);
        win_pop = K_m;
        r_meas = r_meas(:,(K-K_m+1):K);
        if use_odom
            del_X_o = del_X_o(:,:,(K-K_m+1):K-1);
        end
        
        % Re-initialize new W_Q and W_R matrices using most up-to-date 
        % P_check_0 and remaining non-marginalized states
        W_Q_swf_span = (K-K_m+1)*num_xi+1:(K)*num_xi;
        W_R_swf_span = (K-K_m)*length(R)+1:(K)*length(R);
        W_Q = blkdiag(P_check_0\eye(size(P_check_0)), W_Q(W_Q_swf_span,W_Q_swf_span));
        W_R = W_R(W_R_swf_span,W_R_swf_span);
        
        if use_odom
            W_Q_o_swf_span = (K-K_m)*length(Sigma_del_X_o)+1:(K-1)*length(Sigma_del_X_o);
            W_Q_o = W_Q_o(W_Q_o_swf_span,W_Q_o_swf_span);
        end
    end
end

%% Compute final RMSE
RMSE = sqrt(mean(mean(err_batch.^2)))

for ii = 1:length(err_batch)
    RMSE_total(ii) = sqrt(mean(mean(err_batch(:,1:ii).^2)));
    RMSE_yaw(ii) = sqrt(mean(mean(err_batch(3,1:ii).^2)));
    RMSE_pos(ii) = sqrt(mean(mean(err_batch(7:9,1:ii).^2)));
    
    NEES_total(ii) = err_batch(:,ii)' / Sig_batch(:,:,ii) * err_batch(:,ii);
    NEES_yaw(ii) = err_batch(3,ii)' / Sig_batch(3,3,ii) * err_batch(3,ii);
    NEES_pos(ii) = err_batch(7:9,ii)' / Sig_batch(7:9,7:9,ii) * err_batch(7:9,ii);
end

if error_deff == "left"
    RMSE_total_LI = RMSE_total;
    RMSE_yaw_LI = RMSE_yaw;
    RMSE_pos_LI = RMSE_pos;
    NEES_total_LI = NEES_total;
    NEES_yaw_LI = NEES_yaw;
    NEES_pos_LI = NEES_pos;
else
    RMSE_total_RI = RMSE_total;
    RMSE_yaw_RI = RMSE_yaw;
    RMSE_pos_RI = RMSE_pos;
    NEES_total_RI = NEES_total;
    NEES_yaw_RI = NEES_yaw;
    NEES_pos_RI = NEES_pos;
end
        
%% Function to evaluate Jacobians and errors
function [F_1, F_2, G_1] = make_err_Jacobians(X_i, del_X, Del_t, bias_J, b_bar, r_meas, err_deff)
    K = size(X_i, 3);

    if err_deff == "left"
        F_1 = [];
        F_2 = form_F_2_LEFT();       % F_2 has K constant terms so form first one
        G_1 = form_G_1_LEFT();       % G_1 has K constant terms so form first one
    else
        F_2 = form_F_2_RIGHT();       
        G_1 = form_G_1_RIGHT(X_i(:,:,1));
    end
    
    for ii = 1:K-1
        if err_deff == "left"
            F_1(:,:,ii) = form_F_1_LEFT(X_i(:,:,ii+1), X_i(:,:,ii), del_X(:,:,ii), Del_t(ii), bias_J(:,:,:,ii));
            F_2(:,:,ii+1) = form_F_2_LEFT();
            G_1(:,:,ii+1) = form_G_1_LEFT();
        else
            F_1(:,:,ii) = form_F_1_RIGHT(X_i(:,:,ii+1), X_i(:,:,ii), del_X(:,:,ii), Del_t(ii), bias_J(:,:,:,ii));
            F_2(:,:,ii+1) = form_F_2_RIGHT();
            G_1(:,:,ii+1) = form_G_1_RIGHT(X_i(:,:,ii+1));
        end
        
        % Gravity vector
        g_a = [0; 0; -9.81];
        X_bar_1 = X_i(:,:,ii);
        X_bar_2 = X_i(:,:,ii+1);
        F_1complexstep = zeros(15,15);
        F_2complexstep = zeros(15,15);
        F_3complexstep = zeros(3,15);
        h = 1e-5;   % NOTE: Small bug in SO(3) tools, choose h > 1e-6. 
        for lv1 = 1 : 15
            % Complex step F_1
            h_i_c      = zeros(15, 1);
            h_i_c(lv1) = h * 1i;
            if err_deff == "left"
                X_i_c1      = X_bar_1 * se23xR3x2alg.expMap(-h_i_c);
            else
                X_i_c1      = se23xR3x2alg.expMap(-h_i_c) * X_bar_1;
            end
            e_i_c1      = form_e_u_MLG(X_i(:,:,ii+1), X_i_c1, del_X(:,:,ii), Del_t(ii), bias_J(:,:,:,ii), b_bar(:,ii), g_a, err_deff);
            
            F_1complexstep(:, lv1) = (1/h) .* imag(e_i_c1);
            
            % Complex step F_2
            if err_deff == "left"
                X_i_c2      = X_bar_2 * se23xR3x2alg.expMap(-h_i_c);
            else
                X_i_c2      = se23xR3x2alg.expMap(-h_i_c) * X_bar_2;
            end
            e_i_c2      = form_e_u_MLG(X_i_c2, X_i(:,:,ii), del_X(:,:,ii), Del_t(ii), bias_J(:,:,:,ii), b_bar(:,ii), g_a, err_deff);
            F_2complexstep(:, lv1) = (1/h) .* imag(e_i_c2);
            
            % Complex step G_1
            e_i_c3      = form_e_y(X_i_c2, r_meas(:,ii), err_deff);
            F_3complexstep(:, lv1) = (1/h) .* imag(e_i_c3);
        end

        %form_F_1(X_i(:,:,ii+1), X_i(:,:,ii), del_X(:,:,ii), Del_t(ii), bias_J(:,:,:,ii))
        %F_1complexstep
        %norm(F_3complexstep - G_1(:,:,ii+1));
        F_1(:,:,ii) = F_1complexstep;
        F_2(:,:,ii+1) = F_2complexstep;
        %G_1(:,:,ii+1) = F_3complexstep;
    end
end

function F1 = form_F_1_LEFT(X_k, X_kmin, del_X, Del_t, bias_J)
    % Get components
    [C_kmin, v_kmin, r_kmin, b_g_kmin, b_a_kmin] = SE23xR3x2.decompose(X_kmin);
    [C_k, v_k, r_k, b_g_k, b_a_k] = SE23xR3x2.decompose(X_k);
    [del_C, del_v, del_r, del_b_g, del_b_a] = SE23xR3x2.decompose(del_X);
    [del_C_b_g, del_v_b_a, del_v_b_g, del_r_b_a, del_r_b_g] = disassemble_bias_J(bias_J);
    
    J_inv = SO3.computeJLeftInv( SO3.decompose(C_k' * C_kmin * del_C) );
    Z = J_inv * C_k' * C_kmin;
    F1 = [-del_C',                zeros(3,3),   zeros(3,3), -del_C_b_g,   zeros(3,3);
          Z*so3alg.cross(del_v),   -Z,          zeros(3,3), -Z*del_v_b_g, -Z*del_v_b_a;
          Z*so3alg.cross(del_r),   -Z*Del_t,      -Z,       -Z*del_r_b_g,  -Z*del_r_b_a;
          zeros(3,3),             zeros(3,3),   zeros(3,3),  -eye(3),    zeros(3,3);
          zeros(3,3),             zeros(3,3),   zeros(3,3),  zeros(3,3), -eye(3)];
end

function F2 = form_F_2_LEFT()
    F2 = eye(15);
end

function G1 = form_G_1_LEFT()
    G1 = [zeros(3,3), zeros(3,3), eye(3,3), zeros(3,3), zeros(3,3)];
end

function F1 = form_F_1_RIGHT(X_k, X_kmin, del_X, Del_t, bias_J)
    % Get components
    [C_kmin, v_kmin, r_kmin, b_g_kmin, b_a_kmin] = SE23xR3x2.decompose(X_kmin);
    [C_k, v_k, r_k, b_g_k, b_a_k] = SE23xR3x2.decompose(X_k);
    [del_C, del_v, del_r, del_b_g, del_b_a] = SE23xR3x2.decompose(del_X);
    [del_C_b_g, del_v_b_a, del_v_b_g, del_r_b_a, del_r_b_g] = disassemble_bias_J(bias_J);
    
    J_inv = SO3.computeJLeftInv( SO3.decompose(C_kmin * del_C * C_k') );
    Z = J_inv * C_k' * C_kmin;
    v_C = J_inv * so3alg.cross(v_kmin + C_kmin*del_v + v_k);
    v_b_g = J_inv * (C_kmin * del_v_b_g + so3alg.cross(v_k)*C_k*del_C_b_g);
    v_b_a = J_inv * C_kmin * del_v_b_a;
    F1 = [-eye(3),                zeros(3,3),   zeros(3,3), -C_k*del_C_b_g,   zeros(3,3);
          v_C,   -J_inv,          zeros(3,3),                  -v_b_g,         -v_b_a;
          zeros(3,3),   zeros(3,3),      zeros(3,3),       zeros(3,3),  zeros(3,3);
          zeros(3,3),             zeros(3,3),   zeros(3,3),  -eye(3),    zeros(3,3);
          zeros(3,3),             zeros(3,3),   zeros(3,3),  zeros(3,3), -eye(3)];
end

function F2 = form_F_2_RIGHT()
    F2 = eye(15);
end

function G1 = form_G_1_RIGHT(X_k)
    [~, ~, r_k, ~, ~] = SE23xR3x2.decompose(X_k);
    G1 = [-so3alg.cross(r_k), zeros(3,3), eye(3,3), zeros(3,3), zeros(3,3)];
end

function [G_o_1, G_o_2] = make_odom_err_Jacobians(X_i, del_X_o, error_deff)
    K = size(X_i, 3);
    G_o_1 = [];
    G_o_2 = [];
    for ii = 1:K-1
        G_o_1(:,:,ii) = form_G_o_1(X_i(:,:,ii+1), X_i(:,:,ii), del_X_o(:,:,ii));
        G_o_2(:,:,ii) = form_G_o_2(X_i(:,:,ii+1), X_i(:,:,ii), del_X_o(:,:,ii));
        
        % Complex Jacobian eval
        X_bar_1 = X_i(:,:,ii);
        X_bar_2 = X_i(:,:,ii+1);
        G_o_1complexstep = zeros(6,15);
        G_o_2complexstep = zeros(6,15);
        
        h = 1e-5;   % NOTE: Small bug in SO(3) tools, choose h > 1e-6. 
        for lv1 = 1 : 15
            % Complex step F_1
            h_i_c       = zeros(15, 1);
            h_i_c(lv1)  = h * 1i;
            if error_deff == "left"
                X_i_c1      = X_bar_1 * se23xR3x2alg.expMap(-h_i_c);
            else
                X_i_c1      = se23xR3x2alg.expMap(-h_i_c) * X_bar_1;
            end
            e_i_c1      = form_e_o_MLG(X_i(:,:,ii+1), X_i_c1, del_X_o(:,:,ii), error_deff);
            
            G_o_1complexstep(:, lv1) = (1/h) .* imag(e_i_c1);
            
            % Complex step F_2
            if error_deff == "left"
                X_i_c2      = X_bar_2 * se23xR3x2alg.expMap(-h_i_c);
            else
                X_i_c2      = se23xR3x2alg.expMap(-h_i_c) * X_bar_2;
            end
            e_i_c2      = form_e_o_MLG(X_i_c2, X_i(:,:,ii), del_X_o(:,:,ii), error_deff);
            G_o_2complexstep(:, lv1) = (1/h) .* imag(e_i_c2);
        end
        
        G_o_1(:,:,ii) = G_o_1complexstep;
        G_o_2(:,:,ii) = G_o_2complexstep;
        %form_F_1(X_i(:,:,ii+1), X_i(:,:,ii), del_X(:,:,ii), Del_t(ii), bias_J(:,:,:,ii))
        %F_1complexstep
        %F_1(:,:,ii) = F_1complexstep;
    end
end

function Go1 = form_G_o_1(X_k, X_kmin, del_X_o)
    % Get components
    [C_kmin, ~, r_kmin, ~, ~] = SE23xR3x2.decompose(X_kmin);
    [C_k, ~, r_k, ~, ~] = SE23xR3x2.decompose(X_k);
    [del_C_o, del_r_o] = SE3.decompose(del_X_o);
    
    J_inv_phi = SO3.computeJLeftInv( SO3.decompose(del_C_o' * C_kmin' * C_k) );

    Go1 = [J_inv_phi*del_C_o', zeros(3,12);
           -J_inv_phi*del_C_o'*so3alg.cross(del_r_o), zeros(3,3), J_inv_phi*del_C_o', zeros(3,6)];
end

function Go2 = form_G_o_2(X_k, X_kmin, del_X_o)
    [C_kmin, ~, ~, ~, ~] = SE23xR3x2.decompose(X_kmin);
    [C_k, ~, ~, ~, ~] = SE23xR3x2.decompose(X_k);
    [del_C_o, ~] = SE3.decompose(del_X_o);
    J_inv = SO3.computeJRightInv( SO3.decompose(del_C_o' * C_kmin' * C_k) );
    Go2 = [-J_inv*eye(3), zeros(3,12);
           zeros(3,6), -J_inv*eye(3), zeros(3,6)];
end

function H = make_H(F_1, F_2, G_1, G_o_1, G_o_2, use_odom)
    % Construct H matrix
    K = size(F_2,3);
    num_xi = size(F_2,2);
    num_meas = size(G_1,1);
    num_od = size(G_o_1,1);
    H_top = zeros(num_xi*K);
    H_bot = zeros(num_meas*K, num_xi*K);
    
    if use_odom
        H_mid = zeros(num_od*(K-1), num_xi*K);
    end
    
    for ii = 1:K
        % Compose top of H
        H_top((ii-1)*num_xi+1:ii*num_xi, (ii-1)*num_xi+1:ii*num_xi) = F_2(:,:,ii);
        if ii < K
            H_top(ii*num_xi+1:(ii+1)*num_xi, (ii-1)*num_xi+1:ii*num_xi) = F_1(:,:,ii);
        end
        
        if use_odom && ii < K
            % Compose mid of H
            H_mid((ii-1)*(num_od)+1:ii*(num_od), (ii-1)*(num_xi)+1:ii*(num_xi)) = G_o_1(:,:,ii);
            H_mid((ii-1)*(num_od)+1:ii*(num_od), ii*num_xi+1:(ii+1)*num_xi) = G_o_2(:,:,ii);
        end

        % Compose bot of H
        H_bot((ii-1)*num_meas+1:ii*num_meas, (ii-1)*(num_xi)+1:ii*(num_xi)) = G_1(:,:,ii);
    end
    
    if use_odom
        H = [H_top; H_mid; H_bot];
    else
        H = [H_top; H_bot];
    end
end

function e = make_e(X_i, X_check_0, del_X, Del_t, bias_J, b_bar, g_a, r_meas, K, num_xi, del_X_o, use_odom, error_deff)
    % Run trial error to get sizes
    if use_odom
        e_o_t = form_e_o_MLG(X_i(:,:,2), X_i(:,:,1), del_X_o(:,:,1), error_deff);
        num_od = length(e_o_t);
    end
    e_y_t = form_e_y(X_i(:,:,1), r_meas(:,1), error_deff);
    num_meas = length(e_y_t);
    
    % Get initial error
    e_0 = form_e_0_MLG(X_i(:,:,1), X_check_0, error_deff);
    % Initialize process and measurement error
    e_u = zeros([(K-1)*num_xi,1]);
    e_y = zeros([num_meas*K,1]);
    
    if use_odom
        e_o = zeros([num_od*(K-1),1]);
    end
    
    for ii = 1:K
        % Get process error
        if ii > 1
            e_u((ii-2)*num_xi+1:(ii-1)*num_xi) = form_e_u_MLG(X_i(:,:,ii), X_i(:,:,ii-1), del_X(:,:,ii-1), Del_t(ii-1), bias_J(:,:,:,ii-1), b_bar(:,ii-1), g_a, error_deff);
        end
        
        if use_odom && ii > 1
            % Get odometry error
            e_o((ii-2)*num_od+1:(ii-1)*num_od) = form_e_o_MLG(X_i(:,:,ii), X_i(:,:,ii-1), del_X_o(:,:,ii-1), error_deff);
        end
        
        % Get measurement error
        e_y((ii-1)*num_meas+1:ii*num_meas) = form_e_y(X_i(:,:,ii), r_meas(:,ii), error_deff);
    end
    
    if use_odom
        e = [e_0; e_u; e_o; e_y];
    else
        e = [e_0; e_u; e_y];
    end
end

function e_0 = form_e_0_MLG(X_0, X_check_0, err_deff)
    if err_deff == "left"
        e_0 = se23xR3x2alg.vee( SE23xR3x2.logMap( X_0 \ X_check_0 ) );
    else
        e_0 = se23xR3x2alg.vee( SE23xR3x2.logMap( X_check_0 / X_0 ) );
    end
end

function e_u = form_e_u_MLG(X_k, X_kmin, del_X_kmin_k, Del_t_kmin_k, bias_J, b_bar, g_a, err_deff)
    % Update RMIs with new bias
    [~, ~, ~, b_i_g, b_i_a] = SE23xR3x2.decompose(X_kmin);
    b_i = [b_i_g; b_i_a];
    del_X_kmin_k_prop = prop_bias(del_X_kmin_k, bias_J, b_i, b_bar);
    % Get next step predicted through process model
    X_k_check = F(X_kmin, del_X_kmin_k_prop, Del_t_kmin_k, g_a);
    % Get error on full state
    
    if err_deff == "left"
        E_u = X_k \ X_k_check;
    else
        E_u = X_k_check / X_k;
    end
    
    % Get R^n error
    e_u = se23xR3x2alg.vee(SE23xR3x2.logMap(E_u));
end

function e_u = form_e_u_comp(X_k, X_kmin, del_X_kmin_k, Del_t_kmin_k, bias_J, b_bar, g_a, err_deff)
    % Update RMIs with new bias
    [~, ~, ~, b_i_g, b_i_a] = SE23xR3x2.decompose(X_kmin);
    b_i = [b_i_g; b_i_a];
    del_X_kmin_k_prop = prop_bias(del_X_kmin_k, bias_J, b_i, b_bar);
    
    % Get next step predicted through process model
    X_k_check = F(X_kmin, del_X_kmin_k_prop, Del_t_kmin_k, g_a);
    [C_check, v_check, r_check, b_g_check, b_a_check] = SE23xR3x2.decompose(X_k_check);
    
    % Get components of predicted next state
    [C_k, v_k, r_k, b_g_k, b_a_k] = SE23xR3x2.decompose(X_k);
    
    % Get DCM component error
    if err_deff == "left"
        E_phi = C_k' * C_check;
        e_phi = so3alg.vee( SO3.logMap(E_phi) );
        % Get J inverse for other errors
        J_inv = SO3.computeJLeftInv( e_phi );
        % Get velocity component error
        E_v = C_k' * ( v_check - v_k );
        e_v = J_inv * E_v;
        % Get position component error
        E_r = C_k' * ( r_check - r_k );
        e_r = J_inv * E_r;
    else
        E_phi = C_check * C_k';
        e_phi = so3alg.vee( SO3.logMap(E_phi) );
        % Get J inverse for other errors
        J_inv = SO3.computeJLeftInv( e_phi );
        % Get velocity component error
        E_v = v_check - C_check * C_k' * v_k;
        e_v = J_inv * E_v;
        % Get position component error
        E_r = r_check - C_check * C_k' * r_k;
        e_r = J_inv * E_r;
    end

    % Get bias component errors
    e_b_g = b_g_check - b_g_k;
    e_b_a = b_a_check - b_a_k;
    
    % Assemble error terms
    e_u = [e_phi; e_v; e_r; e_b_g; e_b_a];
end

function e_y = form_e_y(X_i, r_meas, err_deff)
    [C_i, ~, r_i, ~, ~] = SE23xR3x2.decompose(X_i);
    if err_deff == "left"
        e_y = C_i' * (r_meas - r_i);
    else
        e_y = r_meas - r_i;
    end
end

function e_o = form_e_o_MLG(X_k, X_kmin, del_X_o, err_deff)
    % Get predicted RMI measurement
    Y_check = G_o(X_k, X_kmin);
    
    if err_deff == "left"
        E_o = del_X_o \ Y_check;
    else
        E_o = Y_check / del_X_o;
    end
    
    e_o = se3alg.vee( SE3.logMap( E_o ) );
end

%% Process model function
function X_k_1 = F(X_k, del_X, Del_t, g_a)
    % Extract components from state
    [C, v, r, b_g, b_a] = SE23xR3x2.decompose(X_k);
    % Make process model matrices
    Gamma = SE23xR3x2.synthesize(eye(3), g_a*Del_t, 0.5*g_a*Del_t^2, zeros(3,1), zeros(3,1));
    Phi = SE23xR3x2.synthesize(C, v, (r+v*Del_t), b_g, b_a);
    X_k_1 = Gamma*Phi*del_X;
end

function X_k_1 = F_comp(X_k, del_X, Del_t, g_a)
    % Extract components from state
    [C_k, v_k, r_k, b_g_k, b_a_k] = SE23xR3x2.decompose(X_k);
    [del_C, del_v, del_r, ~, ~] = SE23xR3x2.decompose(del_X);
    % Update process model
    C_k_1 = C_k * del_C;
    v_k_1 = v_k + Del_t*g_a + C_k*del_v;
    r_k_1 = r_k + Del_t*v_k + 0.5*Del_t^2*g_a + C_k*del_r;
    b_g_k_1 = b_g_k;
    b_a_k_1 = b_a_k;
    
    X_k_1 = SE23xR3x2.synthesize(C_k_1, v_k_1, r_k_1, b_g_k_1, b_a_k_1);
end

function Y_check = G_o(X_k, X_kmin)
    Q = [ eye(3,3) , zeros(3,1), zeros(3,1), zeros(3,5);
         zeros(1,3),    0      ,      1    , zeros(1,5)];
    R = [ eye(3,3) ,  zeros(3,1);
         zeros(1,3),     0      ;
         zeros(1,3),     1      ;
         zeros(5,3),  zeros(5,1)];
    % Get predicted RMI measurement
    Y_check = (Q*X_kmin*R) \ (Q*X_k*R);
end


%% Preintegration
function [del_X, Sigma, Del_t, bias_J, b_bar] = preint_IMU(X_bar, gyr_dat, a_dat, Q_c, dt, error_deff)
    % Preintegrate IMU
    % Get bar state terms
    [C_bar, v_bar, r_bar, b_g_bar, b_a_bar] = SE23xR3x2.decompose(X_bar);
    
    % Save nominal biases
    b_bar = [b_g_bar; b_a_bar];
    
    del_C = eye(3);
    del_v = zeros(3,1);
    del_r = zeros(3,1);
    del_b = zeros(3,1);
    Sigma = zeros(15);  % Hard coded size, need to update value with first Q_d
    Del_t = 0;
    del_C_b_g_2 = zeros(3);
    del_v_b_a = zeros(3);
    del_v_b_g = zeros(3);
    del_r_b_a = zeros(3);
    del_r_b_g = zeros(3);
    % To emulate real problem, loop through IMU measurements as though
    % they're arriving one by one even though we have access to all.
    for ii = 1:size(gyr_dat, 2)
        if ii < size(gyr_dat, 2)
            gyr_dat_ii = (gyr_dat(:,ii) + gyr_dat(:,ii+1))/2;
            a_dat_ii = (a_dat(:,ii) + a_dat(:,ii+1))/2;
        else
            gyr_dat_ii = gyr_dat(:,ii);
            a_dat_ii = a_dat(:,ii);
        end
        % Compute bias Jacobians
        % DCM Jacobian is confusing so just redo every time
        del_C_b_g = zeros(3);
        for kk = 1:ii
            if kk < ii
                gyr_dat_kk = (gyr_dat(:,kk) + gyr_dat(:,kk+1))/2;
            else
                gyr_dat_kk = gyr_dat(:,kk);
            end
            J_r_k = SO3.computeJRight( (gyr_dat_kk -  b_g_bar) * dt );
            del_C_k1_j = compute_del_C(gyr_dat_kk, b_g_bar, dt);
            del_C_b_g = del_C_b_g - del_C_k1_j' * J_r_k * dt;
        end
        %J_r_k = SO3.computeJRight( (gyr_dat(:,ii) -  b_g_bar) * dt );
        %del_C_b_g_2 = del_C_b_g_2 - del_C' * J_r_k * dt;
        %del_C_b_g = del_C_b_g_2;
        del_v_b_a = del_v_b_a - del_C * dt;
        del_v_b_g = del_v_b_g - (del_C * so3alg.cross(a_dat_ii - b_a_bar) * del_C_b_g * dt);
        del_r_b_a = del_r_b_a + (del_v_b_a*dt - 0.5 * del_C * dt^2);
        del_r_b_g = del_r_b_g + (del_v_b_g*dt - 0.5* del_C * so3alg.cross(a_dat_ii - b_a_bar) * del_C_b_g * dt^2);
        
        % Get RMIs
        del_r = del_r + del_v * dt + 0.5*( del_C * (a_dat_ii - b_a_bar) * dt^2 );
        del_v = del_v + del_C * (a_dat_ii - b_a_bar) * dt;
        del_C = del_C*so3alg.expMap(so3alg.cross( (gyr_dat_ii -  b_g_bar) * dt));
        
        % Get uncertainty
        if error_deff == "left"
            [F_d, Q_d] = get_RMI_jacobians_LEFT(gyr_dat_ii, a_dat_ii, b_g_bar, b_a_bar, Q_c, dt);
        else
            [F_d, Q_d] = get_RMI_jacobians_RIGHT(del_C, del_v, del_r, Q_c, dt);
        end
        % Initialize Sigma correctly
        Sigma = F_d * Sigma * F_d' + Q_d;
        
        % Get Delta t
        Del_t = Del_t + dt;
    end
    
    bias_J = assemble_bias_J(del_C_b_g, del_v_b_a, del_v_b_g, del_r_b_a, del_r_b_g);
    
    % Assemble final RMI term
    del_X = SE23xR3x2.synthesize(del_C, del_v, del_r, del_b, del_b);
end

function del_C = compute_del_C(u_g, b_g_bar, dt)
    del_C = eye(3);
    for ii = 1:size(u_g,2)
        del_C = del_C*so3alg.expMap(so3alg.cross( (u_g(:,ii) -  b_g_bar) * dt));
    end
end

function [del_X_o, Sigma] = preint_odom(o_lin_dat, o_ang_dat, Q_o_c, dt, error_deff)
    del_C_o = eye(3);
    del_r_o = zeros(3,1);
    Sigma = zeros(6);  % Hard coded size, need to update value with first Q_d
    % To emulate real problem, loop through odom measurements as though
    % they're arriving one by one even though we have access to all.
    for ii = 1:size(o_lin_dat, 2)
        if ii < size(o_lin_dat, 2)
            o_lin_dat_ii = (o_lin_dat(:,ii) + o_lin_dat(:,ii+1))/2;
            o_ang_dat_ii = (o_ang_dat(:,ii) + o_ang_dat(:,ii+1))/2;
        else
            o_lin_dat_ii = o_lin_dat(:,ii);
            o_ang_dat_ii = o_ang_dat(:,ii);
        end
        
        % Propagate RMIs
        del_r_o = del_r_o + del_C_o * o_lin_dat_ii * dt;
        del_C_o = del_C_o * so3alg.expMap( so3alg.cross( o_ang_dat_ii * dt ) );
        
        % Propagate uncertainty
        if error_deff == "left"
            [F_o_d, Q_o_d] = get_odom_RMI_jacobians_LEFT(o_lin_dat_ii, o_ang_dat_ii, Q_o_c, dt);
        else
            [F_o_d, Q_o_d] = get_odom_RMI_jacobians_RIGHT(del_C_o, del_r_o, Q_o_c, dt);
        end
        % Initialize Sigma correctly
        Sigma = F_o_d * Sigma * F_o_d' + Q_o_d;
    end
    % Assemble final element
    del_X_o = SE3.synthesize(del_C_o, del_r_o);
end

function [F_d, Q_d] = get_RMI_jacobians_LEFT(u_g, u_a, b_g, b_a, Q_c, dt)
    % Make continuous time F_r
    F_r = [-so3alg.cross(u_g - b_g), zeros(3,6), -eye(3,3), zeros(3,3);
           -so3alg.cross(u_a - b_a), -so3alg.cross(u_g - b_g), zeros(3,6), -eye(3,3);
           zeros(3,3), eye(3,3), -so3alg.cross(u_g - b_g), zeros(3,6);
           zeros(6,15)];
    
    % Make continuous time L_r which is constant
    L_r = [eye(3,3), zeros(3,9);
           zeros(3,3), eye(3,3), zeros(3,6);
           zeros(3,12);
           zeros(3,6), -eye(3,3), zeros(3,3);
           zeros(3,9), -eye(3,3)];
       
    [F_d, Q_d] = van_loans(F_r, L_r, Q_c, dt);
end

function [F_d, Q_d] = get_RMI_jacobians_RIGHT(del_C, del_v, del_r, Q_c, dt)
    % Make continuous time F_r
    F_r = [zeros(3,9), -del_C, zeros(3,3);
           zeros(3,9), -so3alg.cross(del_v)*del_C, -del_C;
           zeros(3,3), eye(3,3), zeros(3,3), -so3alg.cross(del_r)*del_C, zeros(3,3);
           zeros(6,15)];
    
    % Make continuous time L_r which is constant
    L_r = [del_C, zeros(3,9);
           so3alg.cross(del_v)*del_C, del_C, zeros(3,6);
           so3alg.cross(del_r)*del_C, zeros(3,9);
           zeros(3,6), -eye(3,3), zeros(3,3);
           zeros(3,9), -eye(3,3)];
       
    [F_d, Q_d] = van_loans(F_r, L_r, Q_c, dt);
end

function [F_o_d, Q_o_d] = get_odom_RMI_jacobians_LEFT(u_v, u_vg, Q_o_c, dt)
    % Make continuous time F_o_r
    F_o_r = [-so3alg.cross(u_vg),      zeros(3,3);
             -so3alg.cross(u_v), -so3alg.cross(u_vg)];
    
    % Make continuous time L_r which is constant
    L_o_r = eye(6);
       
    [F_o_d, Q_o_d] = van_loans(F_o_r, L_o_r, Q_o_c, dt);
end

function [F_o_d, Q_o_d] = get_odom_RMI_jacobians_RIGHT(del_C_o, del_r_o, Q_o_c, dt)
    % Make continuous time F_o_r
    F_o_r = zeros(6);
    
    % Make continuous time L_r which is constant
    L_o_r = [      del_C_o,               zeros(3);
             so3alg.cross(del_r_o),       del_C_o];
       
    [F_o_d, Q_o_d] = van_loans(F_o_r, L_o_r, Q_o_c, dt);
end

%% Helper functions
function t_idx = find_closest_t_idx(timestream, ref_time)
    [d,t_idx] = min(abs(timestream - ref_time));
end

function [A_d, Q_d] = van_loans(A_c, L_c, Q_c, dt)
    sz = size(A_c);
    Psi = [A_c, L_c*Q_c*L_c' zeros(sz), zeros(sz);
           zeros(sz), -A_c', zeros(sz), zeros(sz);
           zeros(sz), zeros(sz), A_c, zeros(sz);
           zeros(sz), zeros(sz), zeros(sz), zeros(sz)];
    Upsilon = expm(Psi*dt);
    Ups_11 = Upsilon(1:sz,1:sz);
    Ups_12 = Upsilon(1:sz,sz+1:2*sz);
    
    A_d = Ups_11;
    Q_d = Ups_12*Ups_11';
end

function del_X = prop_bias(del_X_og, bias_J, b_i, b_bar)
    % Calculate del_b
    del_b = b_i - b_bar;
    % Extract del gyro and accel biases
    del_b_g = del_b(1:3);
    del_b_a = del_b(4:6);
    % Extract individual bias Jacobians and RMI components
    [del_C_b_g, del_v_b_a, del_v_b_g, del_r_b_a, del_r_b_g] = disassemble_bias_J(bias_J);
    [del_C_og, del_v_og, del_r_og, ~, ~] = SE23xR3x2.decompose(del_X_og);
    % Propagate per component
    del_C = del_C_og * so3alg.expMap( so3alg.cross( del_C_b_g*del_b_g) );
    del_v = del_v_og + del_v_b_a*del_b_a + del_v_b_g*del_b_g;
    del_r = del_r_og + del_r_b_a*del_b_a + del_r_b_g*del_b_g;
    % Assemble final 
    del_X = SE23xR3x2.synthesize(del_C, del_v, del_r, zeros(3,1), zeros(3,1));
end

function bias_Jacobians = assemble_bias_J(del_C_b_g, del_v_b_a, del_v_b_g, del_r_b_a, del_r_b_g)
    bias_Jacobians(:,:,1) = del_C_b_g;
    bias_Jacobians(:,:,2) = del_v_b_a;
    bias_Jacobians(:,:,3) = del_v_b_g;
    bias_Jacobians(:,:,4) = del_r_b_a;
    bias_Jacobians(:,:,5) = del_r_b_g;
end

function [del_C_b_g, del_v_b_a, del_v_b_g, del_r_b_a, del_r_b_g] = disassemble_bias_J(bias_Jacobians)
    del_C_b_g = bias_Jacobians(:,:,1);
    del_v_b_a = bias_Jacobians(:,:,2);
    del_v_b_g = bias_Jacobians(:,:,3);
    del_r_b_a = bias_Jacobians(:,:,4);
    del_r_b_g = bias_Jacobians(:,:,5);
end

%% Unit tests
function result = process_error_unit_test()
    % Get MLG-based error
    % Define initial state X_kmin
    X_kmin = SE23R3R3_synthesize(zeros(3,1), zeros(3,1), zeros(3,1), zeros(3,1), zeros(3,1));
    
    % Gravity vector
    g_a = [0; 0; -9.81];
    
    % Gen sensor data
    sen_num = 100;
    dt = 0.01;
    
    sigma_accel = 0.001;      % [m / s^2]
    sigma_gyro = 0.005;      % [rad / s]
    sigma_accel_bias = 0.01; % [m / s^3]
    sigma_gyro_bias = 0.01;  % [rad / s^2]

    Q_r = diag([sigma_gyro^2, sigma_gyro^2, sigma_gyro^2,...
            sigma_accel^2, sigma_accel^2, sigma_accel^2,...
            sigma_gyro_bias^2, sigma_gyro_bias^2, sigma_gyro_bias^2,...
            sigma_accel_bias^2, sigma_accel_bias^2, sigma_accel_bias^2]); % Continuous time PSD
    
    % Gen gt data
    gyr_dat_gt = zeros(3,sen_num) + [0; 0; 1];
    a_dat_gt = zeros(3,sen_num) - g_a;
    
    % Gen noisy data
    gyr_dat = gyr_dat_gt + sigma_gyro.*randn(size(gyr_dat_gt));
    a_dat = a_dat_gt + sigma_accel.*randn(size(a_dat_gt));
    
    % Get gt RMIs
    [del_X_kmin_k_gt, ~, Del_t_kmin_k_gt, ~] = preint_IMU(X_kmin, gyr_dat_gt, a_dat_gt, Q_r, dt);
    
    % Get noisy RMIs
    [del_X_kmin_k, ~, Del_t_kmin_k, ~] = preint_IMU(X_kmin, gyr_dat, a_dat, Q_r, dt);
    
    % Gen gt next state
    X_k_gt = F(X_kmin, del_X_kmin_k_gt, Del_t_kmin_k_gt, g_a);
    
    % Gen random next state
    X_k = SE23xR3x2.synthesize(zeros(3,1), zeros(3,1) + [2;0;0], zeros(3,1), zeros(3,1), zeros(3,1));
    
    X_k_use = X_k;
    
    % Form MLG errors
    e_u_MLG_gt = form_e_u_MLG(X_k_use, X_kmin, del_X_kmin_k_gt, Del_t_kmin_k_gt, g_a);
    e_u_MLG = form_e_u_MLG(X_k_use, X_kmin, del_X_kmin_k, Del_t_kmin_k, g_a);
    
    % Form component-wise errors
    e_u_comp_gt = form_e_u_comp(X_k_use, X_kmin, del_X_kmin_k_gt, Del_t_kmin_k_gt, g_a);
    e_u_comp = form_e_u_comp(X_k_use, X_kmin, del_X_kmin_k, Del_t_kmin_k, g_a);
    
    if norm(e_u_MLG_gt - e_u_comp_gt) < 10^-14 && norm(e_u_MLG - e_u_comp) < 10^-14
        result = 1;
    else
        result = 0;
    end
end