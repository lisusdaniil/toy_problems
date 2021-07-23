% AUTHOR:       Daniil Lisus (260669654)
% TOPIC:        Extended Kalman Filter
% MODULE:       Extended Kalman Filter
% DESCRIPTION:  Filters the provided noisy data to provide an estimate for
%               the position and velocity of a cart in front of a wall. The
%               velocity measurements and the process model are still
%               linear, however the range measurement and the measurement
%               model are non-linear and thus need to be modified for the
%               EKF. Also returns the 3 sigma bound for the position and 
%               innovation term.

% Choose which data to use
if plot_ode
    r_dat = r_ode_cor;
    v_dat = v_ode_cor;
else
    r_dat = r_sin_cor;
    v_dat = v_sin_cor;
end

% Define static matricies
L = 1;
M = 1;
Q = 0.1*vel_var^2;
R = 5*range_var^2;

% Define continuous-time state space system
A_c = 0;
B_c = 1;

% Define discrete-time state space system
% Q here is technically continous Q, but we use cov Q as dummy variable
dt = T(2)-T(1);
xi = [A_c L*Q*L' 0      0;
      0   -A_c'  0      0;
      0     0    A_c  B_c;
      0     0    0      0];

% Since our xi is a nillpotent matrix of degree 2, can compute exponential
% exactly
upsilon = eye(size(xi)) + xi*dt;
A_d = upsilon(1,1);
B_d = upsilon(3,4);

% Define initial values
P_hat = R;
r_norm(1) = 3*sqrt(P_hat);

tick = 0;
i_inov = 1;
for ii=2:length(T)
    tick = tick + 1;
    
    % Prediction step
    x_check = A_d*x_hat(ii-1) + B_d*v_dat(ii-1);
    P_check = A_d*P_hat*A_d' + L*Q*L';
    P_check = 0.5*(P_check + P_check');
    
    if (same_freq == 1 || tick == 10) && ~(same_freq == 2)
        tick = 0;
        
        % Correction step
        C_k = x_check/sqrt(l^2+x_check^2);  % Linearized measurement model about x_check
        S_k = C_k*P_check*C_k' + M*R*M';
        S_k = 0.5*(S_k + S_k');              % Enforce symmetry
        K_k = P_check*C_k'/S_k;             % Kalman gain
        
        % Covariance update
        P_hat = (1 - K_k*C_k)*P_check*(1 - K_k*C_k)' + K_k*M*R*M'*K_k'; 
        rho_k(i_inov) = r_dat(ii) - sqrt(x_check^2 + l^2);    % innovation term
        x_hat(ii) = x_check + K_k*rho_k(i_inov);              % final update
        r_used(i_inov) = sqrt(r_dat(ii) - l^2);
        
        inov_norm(i_inov) = 3*sqrt(S_k);
        T_inov(i_inov) = T(ii);
        i_inov = i_inov + 1;
    else
        x_hat(ii) = x_check;
        P_hat = P_check;
    end
    
    P_hat = 0.5*(P_hat + P_hat');       % Enforce symmetry

    % Calculate 3 sigma norms
    r_norm(ii) = 3*sqrt(P_hat);
end

x_KF = x_hat';
