% AUTHOR:       Daniil Lisus (260669654)
% TOPIC:        Kalman Filter
% MODULE:       Runfile
% DESCRIPTION:  Runfile for the entire Kalman Filter exercise. Generates
%               position and velocity data for a "car" in front of a wall,
%               then corrupt the data and use a Kalman Filter to estimate
%               the position and velocity of the car

%% -----------------------------------------------------------------------
%   Setup
%-------------------------------------------------------------------------

% House cleaning
clear;
clc;
close all;

% Time discreitzation
T = 0:0.05:10;

% Noise generating variance
pos_var = 0.1;
vel_var = 0.01;

% Initial guesses
x_hat(1) = 4.0;

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
KF_Gen_Corrupt_Data

%% -----------------------------------------------------------------------
%   Apply Kalman Filter
%-------------------------------------------------------------------------
KF_Filter

%% -----------------------------------------------------------------------
%   Plot
%-------------------------------------------------------------------------
KF_plots

