% AUTHOR:       Daniil Lisus (260669654)
% TOPIC:        MLG Batch Estimation
% MODULE:       Generate and Corrupt Data
% DESCRIPTION:  Generate a robot 3D trajectory based on provided
%               velocities.

t_train = T;
ideal_data = false;
r_uc_b = [0;0;0];

% Initial conditions for robot
r0_a1 = [0.0; 0.0; 0.0];                                % [m, m , m]          
theta0_a1 = [0.0; 0.0; 0.0];                           % [rad, rad , rad]
omega_a1 = [2.0+0.0*t_train; 6.0+0.0*t_train; 7.0+0.0*t_train];     % [rad/s, rad/s, rad/s]
v_a1 = [8.0+0.0*t_train; 2.0+0.0*t_train; 1.0+0.0*t_train];         % [m/s, m/s, m/s]

test_path = robot_path("Robot","None",r_uc_b,ideal_data,t_train,r0_a1,theta0_a1,omega_a1,v_a1);