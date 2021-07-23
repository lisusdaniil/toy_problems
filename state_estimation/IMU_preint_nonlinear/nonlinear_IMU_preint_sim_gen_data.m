% AUTHOR:       Daniil Lisus (260669654)
% TOPIC:        Nonlinear IMU Preintegration Simulated Data
% MODULE:       Generate and Corrupt Data
% DESCRIPTION:  Generates position and velocity data using a simple 
%               sinusoidal function and by integrating a damped spring-mass 
%               system using ode45. Translate the position data to a range 
%               measurement data assuming that the cart tracks a constant
%               point on the wall l meters above the axis of movement.
%               Then adds noise to both the range and velocity generated.

% Generate position/velocity data from simple sinusoidal function
[x_sin,v_sin] = x_v_sin_data(T);

% Generate position/velocity data from spring-mass model
m = 1;                          % (kg) unit mass
k = 1;                          % (N/m) spring coefficient of motion
u = 0;                          % Change for external force control
b = 0.5;                        % Change for damping control
A = [0 1; -k/m -b/m];           % Undamped spring-mass system
B = [0; 1];
x0 = [1; 0];                    % Oscillate between -1 and 1

options = odeset('RelTol',1e-10,'AbsTol',1e-10);
[t_ode,y_ode] = ode45(@(t,y) func(t,y,u,A,B), T, x0, options);

x_ode = 2 + y_ode(:,1);         % Shift oscillation to be about x=2 m
v_ode = y_ode(:,2);

% Change x axis distance to the range measurement of a point l meters up 
r_sin = sqrt(x_sin.^2 + ell.^2);
r_ode = sqrt(x_ode.^2 + ell.^2);

% Corrupt Data
w_r = sqrt(range_var).*randn(size(T));
w_v = sqrt(vel_var).*randn(size(T));
r_sin_cor = r_sin + w_r';
r_ode_cor = r_ode + w_r';
v_sin_cor = v_sin + w_v';
v_ode_cor = v_ode + w_v';

% Change noisy range measurements back to x distance for plotting
x_sin_cor = sqrt(abs(r_sin_cor.^2 - ell.^2));
x_ode_cor = sqrt(abs(r_ode_cor.^2 - ell.^2));


%-------------------------------------------------------------------------
%   Functions
%-------------------------------------------------------------------------

function [x_sin, v_sin] = x_v_sin_data(T)
    x_sin = 2+sin(T+pi/2)'; % "Wall" is at x=0, oscillate between 1 and 3
    v_sin = cos(T+pi/2)';   % Velocity for this model is just derivative
end

function x_dot = func(t,x,u,A,B)
    x_dot = A*x+B*u;
end
