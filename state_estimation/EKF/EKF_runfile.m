% AUTHOR:       Daniil Lisus (260669654)
% TOPIC:        Extended Kalman Filter
% MODULE:       Runfile
% DESCRIPTION:  Runfile for the entire Extended Kalman Filter exercise. 
%               Generates range and velocity data for a cart in front of a 
%               wall, then corrupt the data and use an Extended Kalman 
%               Filter to estimate the position and velocity of the car.
%               Note that the range measurements are assumed to be taken
%               from a constant point on the wall, above the axis of
%               movement and thus the relation between position and range
%               measurement is non-linear.

%% -----------------------------------------------------------------------
%   Setup
%-------------------------------------------------------------------------

% House cleaning
clear;
clc;
close all;

% Time discreitzation
T = 0:0.05:5.0;

% Noise generating variance
range_var = 0.1;
vel_var = 0.01;

% Initial guesses
x_hat(1) = 3;

% 1 - Generate position/velocity data from spring-mass model using ode45
% 0 - Generate position/velocity data from simple sinusoidal function
plot_ode = 0;

% 2 - Run prediction only
% 1 - Measurements come in at same frequency
% 0 - Range is updated at 1/10th frequency of odometry
same_freq = 1;

%% -----------------------------------------------------------------------
%   Generate Data
%-------------------------------------------------------------------------
EKF_Gen_Corrupt_Data

%% -----------------------------------------------------------------------
%   Apply Kalman Filter
%-------------------------------------------------------------------------
EKF_Filter

%% -----------------------------------------------------------------------
%   Plot
%-------------------------------------------------------------------------
EKF_plots

