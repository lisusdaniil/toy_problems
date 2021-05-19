% AUTHOR:       Daniil Lisus (260669654)
% TOPIC:        Linear Batch Estimation
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
T = 0:0.1:10;

% Noise generating variance
pos_var = 0.1;
vel_var = 0.01;

% Initial guesses
x_0 = 2.0;
P_0 = 1.0^2;

% 1 - Generate position/velocity data from spring-mass model using ode45
% 0 - Generate position/velocity data from simple sinusoidal function
plot_ode = 0;

% 2 - Run prediction only
% 1 - Measurements come in at same frequency
% 0 - Range is updated at 1/10th frequency of odometry
same_freq = 0;

%% -----------------------------------------------------------------------
%   Generate Data
%-------------------------------------------------------------------------
lin_batch_gen_corrupt_data

%% -----------------------------------------------------------------------
%   Apply Kalman Filter
%-------------------------------------------------------------------------
lin_batch_filter

%% -----------------------------------------------------------------------
%   Plot
%-------------------------------------------------------------------------
lin_batch_plots

