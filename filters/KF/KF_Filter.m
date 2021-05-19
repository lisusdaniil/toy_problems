% AUTHOR:       Daniil Lisus (260669654)
% TOPIC:        Kalman Filter
% MODULE:       Kalman Filter
% DESCRIPTION:  Filters the provided noisy data to provide an estimate for
%               the position and velocity of a cart in front of a wall.
%               Also returns the 3 sigma bound for the position, velocity
%               and innovation term.

% Choose which data to use
if plot_ode
    r_dat = r_ode_cor;
    v_dat = v_ode_cor;
else
    r_dat = r_sin_cor;
    v_dat = v_sin_cor;
end

% Define static matricies
C = 1;
L = 1;
M = 1;
Q = vel_var^2;
R = 10*pos_var^2;

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
x0 = r_dat(1);                  % First measured values
P_hat = R;
r_norm(1) = 3*sqrt(P_hat);

tick = 0;
i_inov = 1;
for ii=2:length(T)
    tick = tick + 1;
    
    % Prediction step
    x_check = A_d*x_hat(ii-1) + B_d*v_dat(ii-1);
    P_check = A_d*P_hat*A_d' + L*Q*L';
    
    if (same_freq == 1 || tick == 10) && ~(same_freq == 2)
        tick = 0;
        % Correction step
        S_k = C*P_check*C' + M*R*M';
        K_k = P_check*C'/S_k;                       % Kalman gain
        rho_k(i_inov) = r_dat(ii) - C*x_check;        % innovation term
        r_used(i_inov) = r_dat(ii);
        x_hat(ii) = x_check + K_k*rho_k(i_inov);      % final update
        
        % Covariance update
        P_hat = (1 - K_k*C)*P_check*(1 - K_k*C)' + K_k*M*R*M'*K_k';
        P_hat = 0.5*(P_hat + P_hat');       % Enforce symmetry
        
        inov_norm(i_inov) = 3*sqrt(S_k);
        T_inov(i_inov) = T(ii);
        i_inov = i_inov + 1;
    else
        x_hat(ii) = x_check;
        P_hat = P_check;
    end

    % Calculate 3 sigma norms
    r_norm(ii) = 3*sqrt(P_hat);
end

r_KF = x_hat';
