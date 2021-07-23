% AUTHOR:       Daniil Lisus (260669654)
% TOPIC:        Nonlinear IMU Preintegration Simulated Data
% MODULE:       Runfile
% DESCRIPTION:  

%% -----------------------------------------------------------------------
%   Setup
%-------------------------------------------------------------------------

% House cleaning
clear;
clc;
%close all;

% Time discreitzation
dt = 0.01;
T = 0:dt:20.0;

% Noise generating variance
range_var = 0.1^2;
vel_var = 0.5^2;

% Define position up the wall of point being measured
ell = 2.0;     % m

% Initial guesses
x_check_0 = 2.0;
P_check_0 = 1.0^2;

% 1 - Generate position/velocity data from spring-mass model using ode45
% 0 - Generate position/velocity data from simple sinusoidal function
plot_ode = 0;

% Downsample exteroceptive frequency
downsample_freq = 10;

% Control iteration params
max_iter = 10000;
tol_step = 10^-7;   % Step tolerance criteria
tol_grad = 10^-7;   % Gradient tolerance criteria (for LM)
solver = "LM";      % Choose Gauss Newton ("GN") or Levenberg–Marquardt ("LM") 

% Control IMU preintegration
IMU_preint = 0;

%% -----------------------------------------------------------------------
%   Generate Data
%-------------------------------------------------------------------------
nonlinear_IMU_preint_sim_gen_data

%% -----------------------------------------------------------------------
%   Apply Batch Filter
%-------------------------------------------------------------------------
if solver == "GN"
    nonlinear_IMU_preint_batch_filter_GN
elseif solver == "LM"
    nonlinear_IMU_preint_batch_filter_LM
end

%% -----------------------------------------------------------------------
%   Plot
%-------------------------------------------------------------------------
nonlinear_IMU_preint_batch_plots

