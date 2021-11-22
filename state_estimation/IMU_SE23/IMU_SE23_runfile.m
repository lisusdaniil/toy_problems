% AUTHOR:       Daniil Lisus (260669654)
% TOPIC:        SE2(3) with bias IMU Preintegration Simulated Data
% MODULE:       Runfile
% DESCRIPTION:  

%% -----------------------------------------------------------------------
%   Setup
%-------------------------------------------------------------------------
% House cleaning
clear;
clc;
close all;
%%
% Add paths
addpath(genpath('S:/School/Grad/packages/MATLAB_Packages/trajectory_gen/sim_generator'))
addpath(genpath('S:/School/Grad\Packages/MATLAB_Packages/decar_mocap_tools'))

% GT trajectory speed 
traj.speed = 1.0;

% Initial guess
% Initial states
rot_0 = zeros(3,1);
vel_0 = [traj.speed;0; 0];
pos_0 = zeros(3,1);
b_g_0 = 0.00*ones(3,1);
b_a_0 = 0.00*ones(3,1);
% Uncertainty on initial states
% P_rot_0 = (pi/4)^2*ones(1,3);
% P_vel_0 = (1.0)^2*ones(1,3);
% P_pos_0 = (1.0)^2*ones(1,3);
% P_b_g_0 = (0.1)^2*ones(1,3);
% P_b_a_0 = (0.1)^2*ones(1,3);
P_rot_0 = (0.00001)^2*ones(1,3);
P_vel_0 = (0.00001)^2*ones(1,3);
P_pos_0 = (0.00001)^2*ones(1,3);
P_b_g_0 = (0.00001)^2*ones(1,3);
P_b_a_0 = (0.00001)^2*ones(1,3);

% Assemble full element
X_check_0_init = SE23xR3x2.synthesize(rot_0, vel_0, pos_0, b_g_0, b_a_0);
% Assemble full uncertainty
P_check_0_init = diag([P_rot_0, P_vel_0, P_pos_0, P_b_g_0, P_b_a_0]);

% Control iteration params
max_iter = 200;    % Maximum number of iterations for the algorithm
LM_iter = 5;       % Number of iterations to do before turning on LM
tol_step = 10^-7;   % Step tolerance criteria
tol_grad = 10^-7;   % Gradient tolerance criteria (for LM)
%solver = "LM";     % Choose Gauss Newton ("GN") or Levenberg–Marquardt ("LM") 

% Control IMU preintegration
IMU_preint = 1;

% Control sliding window parameters
use_marg = 0;
win_size = 10;
marg_size = 25;

% Control measurements
use_odom = 0;

verbose = true;

error_deff = "right";

%% -----------------------------------------------------------------------
%   Generate Data
%-------------------------------------------------------------------------
config_data_generator
main_data_generator
% gt_SE23R3R3 contains ground truth states
% sensorData contains noisy sensor measurements
%% -----------------------------------------------------------------------
%   Apply Batch Filter
%-------------------------------------------------------------------------
IMU_SE23_preint_batch_filter_LM

%% -----------------------------------------------------------------------
%   Plot
%-------------------------------------------------------------------------
IMU_SE23_preint_batch_plots

