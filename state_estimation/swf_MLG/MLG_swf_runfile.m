% AUTHOR:       Daniil Lisus (260669654)
% TOPIC:        MLG Batch Estimation
% MODULE:       Runfile
% DESCRIPTION:  

%% -----------------------------------------------------------------------
%   Setup
%-------------------------------------------------------------------------

% House cleaning
clear;
clc;
close all;
parameters;
% Change this path to be to the matrix lie tools
addpath(genpath('S:\School\Grad\Packages\MATLAB_Packages\decar_mocap_tools'))

% Time discreitzation
dt = 0.1;
T = 0:dt:40.0;

% Noise generating variance
% Go to Parameters to change noise values

% Initial guess
X_check_0 = SE3.synthesize(zeros(3,1), zeros(3,1));
P_check_0 = diag([(pi/2)^2,(pi/2)^2,(pi/2)^2,1.0^2,1.0^2,1.0^2]);

% Control iteration params
max_iter = 5000;    % Maximum number of iterations for the algorithm
LM_iter = 20;       % Number of iterations to do before turning on LM
tol_step = 10^-12;   % Step tolerance criteria
tol_grad = 10^-12;   % Gradient tolerance criteria (for LM)
solver = "LM";      % Choose Gauss Newton ("GN") or Levenberg–Marquardt ("LM") 

% Control sliding window parameters
win_size = 50;
marg_size = 45;

% Correct T so that we dont have non solved states
T = T(1:(size(T,2)-mod(size(T,2), win_size)));

%% -----------------------------------------------------------------------
%   Generate Data
%-------------------------------------------------------------------------
MLG_swf_gen_corrupt_data

%% -----------------------------------------------------------------------
%   Apply Kalman Filter
%-------------------------------------------------------------------------
if solver == "GN"
    MLG_swf_filter_GN
elseif solver == "LM"
    MLG_swf_filter_LM
end

%% -----------------------------------------------------------------------
%   Plot
%-------------------------------------------------------------------------
MLG_swf_plots

