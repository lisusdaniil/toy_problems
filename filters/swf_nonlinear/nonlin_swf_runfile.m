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
range_var = 0.3^2;
vel_var = 0.05^2;

% Define position up the wall of point being measured
ell = 2.0;     % m

% Initial guesses
x_check_0 = 2.0;
P_check_0 = 2.0^2;

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

% Control sliding window parameters
win_size = 10;
marg_size = 5;

% Correct T so that we dont have non solved states
T = T(1:(size(T,2)-mod(size(T,2), win_size)));

%% -----------------------------------------------------------------------
%   Generate Data
%-------------------------------------------------------------------------
nonlin_swf_gen_corrupt_data

%% -----------------------------------------------------------------------
%   Apply Kalman Filter
%-------------------------------------------------------------------------
if solver == "GN"
    nonlin_swf_filter_GN
elseif solver == "LM"
    nonlin_swf_filter_LM
end
%% -----------------------------------------------------------------------
%   Plot
%-------------------------------------------------------------------------
nonlin_swf_plots

