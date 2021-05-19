% AUTHOR:       Daniil Lisus (260669654)
% TOPIC:        MLG Batch Estimation
% MODULE:       Parameters
% DESCRIPTION:  

% NOISE
GPS_x_std   = 0.1;           % m^2
GPS_y_std   = 0.1;           % m^2
GPS_z_std   = 0.1;           % m^2
gyro_x_std  = 0.01;               % rad/s
gyro_y_std  = 0.01;               % rad/s
gyro_z_std  = 0.01;               % rad/s
vel_x_std   = 0.01;              % m/s
vel_y_std   = 0.01;              % m/s
vel_z_std   = 0.01;              % m/s

gyro_x_bias = 0.0;                   % rad/s
gyro_y_bias = 0.0;                   % rad/s
gyro_z_bias = 0.0;                   % rad/s
gyro_x_walk = 0.0;                  % rad/s^2
gyro_y_walk = 0.0;                  % rad/s^2
gyro_z_walk = 0.0;                  % rad/s^2