% AUTHOR:       Daniil Lisus (260669654)
% TOPIC:        Linear Batch Estimation
% MODULE:       Linear Batch Estimator
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
B = T(2)-T(1);
C = 1;
L = 1;
M = 1;
Q = vel_var;
R = 10*pos_var;

% Get num of states
K = length(r_dat);

% Construct H matrix
H_top = eye(K*length(A));
for ii = 1:K-1
    H_top(ii*length(A)+1:(ii+1)*length(A), 1+(ii-1)*length(A):(ii)*length(A)) = -A;
end
H_bot = kron(eye(K*length(A)),-C);

H = [H_top; H_bot];

% Construct z matrix
z = [];
z(1:length(A),1) = x_0;

for ii = 2:K
   z(ii:ii+length(A)-1,1) = B*v_dat(ii-1);
   z(K+ii-1:K+ii+length(A)-1,1) = -r_dat(ii-1);
end

% Construct W matrix
Q_copy = repmat({Q\eye(size(Q))}, 1, K-1);
R_copy = repmat({R\eye(size(R))}, 1, K);
W = blkdiag(P_0\eye(size(P_0)), Q_copy{:}, R_copy{:});

% Evaluate x_hat
[m, n] = size(H);
Sigma = (H'*W*H)\eye(n);
r_batch = Sigma*H'*W*z;
