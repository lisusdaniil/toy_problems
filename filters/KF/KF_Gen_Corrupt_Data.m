% AUTHOR:       Daniil Lisus (260669654)
% TOPIC:        Kalman Filter
% MODULE:       Generate and Corrupt Data
% DESCRIPTION:  Generates position and velocity data using a simple 
%               sinusoidal function and by integrating a damped spring-mass 
%               system using ode45. Then adds noise to both the position
%               and velocity generated.

% Generate position/velocity data from simple sinusoidal function
[r_sin,v_sin] = r_v_sin_data(T);

% Generate position/velocity data from spring-mass model
m = 1;                          % (kg) unit mass
k = 1;                          % (N/m) spring coefficient of motion
u = 0;                          % Change for external force control
b = 0.2;                          % Change for damping control
A = [0 1; -k/m -b/m];           % Undamped spring-mass system
B = [0; 1];
x0 = [1; 0];                    % Oscillate between -1 and 1

options = odeset('RelTol',1e-10,'AbsTol',1e-10);
[t_ode,y_ode] = ode45(@(t,y) func(t,y,u,A,B), T, x0, options);

r_ode = 2 + y_ode(:,1);         % Shift oscillation to be about r=2 m
v_ode = y_ode(:,2);

%   Corrupt Data
w_r = sqrt(pos_var).*randn(size(T));
w_v = sqrt(vel_var).*randn(size(T));
r_sin_cor = r_sin + w_r';
r_ode_cor = r_ode + w_r';
v_sin_cor = v_sin + w_v';
v_ode_cor = v_ode + w_v';


%-------------------------------------------------------------------------
%   Functions
%-------------------------------------------------------------------------

function [r_sin, v_sin] = r_v_sin_data(T)
    r_sin = 2+sin(T+pi/2)'; % "Wall" is at r=0, oscillate between 1 and 3
    v_sin = cos(T+pi/2)';   % Velocity for this model is just derivative
end

function x_dot = func(t,x,u,A,B)
    x_dot = A*x+B*u;
end
