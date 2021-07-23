% AUTHOR:       Daniil Lisus (260669654)
% TOPIC:        Nonlinear Batch Estimation
% MODULE:       Runfile
% DESCRIPTION:  

%% -----------------------------------------------------------------------
%   Setup
%-------------------------------------------------------------------------

% House cleaning
clear;
clc;
close all;

% Time discreitzation
dt = 0.1;
T = 0:dt:10.0;

% Noise generating variance
range_var = 0.1;
vel_var = 0.01;

% Define position up the wall of point being measured
ell = 2.0;     % m

% Initial guesses
x_check_0 = 2.0;
P_check_0 = 1.0^2;

% 1 - Generate position/velocity data from spring-mass model using ode45
% 0 - Generate position/velocity data from simple sinusoidal function
plot_ode = 0;

% 2 - Run prediction only
% 1 - Measurements come in at same frequency
% 0 - Range is updated at 1/10th frequency of odometry
same_freq = 0;

% Control iteration params
max_iter = 10000;
tol_step = 10^-7;   % Step tolerance criteria
tol_grad = 10^-7;   % Gradient tolerance criteria (for LM)
solver = "LM";      % Choose Gauss Newton ("GN") or Levenberg–Marquardt ("LM") 

%% -----------------------------------------------------------------------
%   Generate Data
%-------------------------------------------------------------------------
nonlin_batch_gen_corrupt_data

%% -----------------------------------------------------------------------
%   Apply Batch Filter
%-------------------------------------------------------------------------
if solver == "GN"
    nonlin_batch_filter_GN
elseif solver == "LM"
    nonlin_batch_filter_LM
end
%% -----------------------------------------------------------------------
%   Plot
%-------------------------------------------------------------------------
nonlin_batch_plots

