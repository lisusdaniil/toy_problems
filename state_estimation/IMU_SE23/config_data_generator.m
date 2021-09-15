%% Trajectory config
traj.lla0.lat0 = 45.262071;     % [deg]
traj.lla0.lon0 = -81.657709;    % [deg]
traj.lla0.alt0 = 0;
% traj.waypoints = [0, 0, 0; 4 4 0; 6 4 0; 6 2 0; 3 2 0; -1 5 0].';


traj.waypoints = [ 0 0 0 ;
                   10, 0, 0;
                   10, 5, 0;
                   0, 5, 0;
                   0, 10, 0;
                   5, 10, 0;
                   5, -2, 0].';
traj.radius = 1.0;
%traj.speed = 5.0;



%% Sensor data config
% North-east-down geodetic frame, or east-north-up
% geodetic_frame_convention = ['ned'] ['enu']
geodetic_frame_convention = 'enu';

% Use WGS84 gravity model, or keep gravity at 9.81.
use_WGS84 = false;

% Treat local geodetic frame as inertial, or include rotation effects for
% IMU measurements. 
include_earth_rotation = false;

% Gravity vector
g_a = [0; 0; -9.81];

%% Sensor config
num_samples = 1;

% Position sensor
f_position = 5;          % [Hz]
sigma_position = 0.01;    % [m]

% IMU
f_imu = 100;             % [Hz]
sigma_accel = 0.001;      % [m / s^2]
sigma_gyro = 0.005;      % [rad / s]
sigma_accel_bias = 0.01; % [m / s^3]
sigma_gyro_bias = 0.01;  % [rad / s^2]
initial_accel_bias = 0.00*ones(3, 1);
initial_gyro_bias = 0.00*ones(3, 1);

% Odometry sensor
f_odom = 100;            % [Hz]
sigma_odom_lin = 0.1;    % [m / s]
sigma_odom_ang = 0.1;    % [rad / s]