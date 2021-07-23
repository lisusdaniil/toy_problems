% Set up filenames
anchor_mocap_file = 'data\_slash_vrpn_client_node_slash_uwb_anchor_linear_slash_pose.csv';
husky_mocap_file = 'data\_slash_vrpn_client_node_slash_Husky_20210607_slash_pose.csv';
odometry_file = 'data\_slash_odometry_slash_filtered.csv';
pozyx_file = 'data\_slash_pozyx_slash_0x6a7b_slash_range.csv';

% Load in real data
anchor_mocap_data = readtable(anchor_mocap_file, 'PreserveVariableNames', true);
husky_mocap_data = readtable(husky_mocap_file, 'PreserveVariableNames', true);
odometry_data = readtable(odometry_file, 'Delimiter', [','], 'PreserveVariableNames', true);
pozyx_data = readtable(pozyx_file, 'PreserveVariableNames', true);

%% Process relevant data

% Set mocap to be t_0
t_0 = husky_mocap_data.rosbagTimestamp(1);

% Get anchor position
anchor_pos = mean(anchor_mocap_data.x); % [m]

% Get ground truth robot state
gt_t = (husky_mocap_data.rosbagTimestamp-t_0)*10^-9; % [s]
% Fit spline to groundtruth
x_gt_c = fit(gt_t, husky_mocap_data.x, 'smoothingspline','SmoothingParam',0.99999999);

%{
figure
plot(x_gt_c, gt_t, husky_mocap_data.x)
RMSE_ = sqrt(mean((x_gt_c(gt_t) - husky_mocap_data.x).^2))
%}
%%
% Get interoceptive measurements
odom_t = (odometry_data.rosbagTimestamp-t_0)*10^-9; % [s]

% Get covariance
odom_cov = [];
for ii = 1:length(odometry_data.covariance_1)
    tmp_cell = odometry_data.covariance_1(ii);
    tmp_mat = cell2mat(tmp_cell);
    tmp_split = strsplit(tmp_mat,{',','['});
    odom_cov(ii) = str2double(cell2mat(tmp_split(2)));
end
odom_cov = odom_cov';

odom_meas = [odom_t, odometry_data.x_2, odom_cov];            % [s, m/s]
odom_meas = odom_meas(odom_meas(:,1)>0,:);          % Delete t < 0 measurements

% Get exteroceptive measurements
range_t = (pozyx_data.rosbagTimestamp-t_0)*10^-9;   % [s]
range_meas = [range_t, pozyx_data.range/1000];      % [s, m]
range_meas = range_meas(range_meas(:,1)>0,:);       % Delete t < 0 measurements



