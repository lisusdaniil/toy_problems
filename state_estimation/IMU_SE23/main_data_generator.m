%% Instantiate sensor objects
myPositionSensor = PositionSensor(f_position, sigma_position);
myIMU = InertialMeasurementUnit(f_imu, sigma_accel, sigma_gyro, ...
    sigma_accel_bias, sigma_gyro_bias, initial_accel_bias, ...
    initial_gyro_bias);
myOdomSensor = OdometrySensor(f_odom, sigma_odom_lin, sigma_odom_ang);

%% Generate the vehicle states
states_position = generateStates(traj, 1 / myPositionSensor.getFrequency(), ...
    geodetic_frame_convention, use_WGS84);
states_imu      = generateStates(traj, 1 / myIMU.getFrequency(), ...
    geodetic_frame_convention, use_WGS84);
states_odom      = generateStates(traj, 1 / myOdomSensor.getFrequency(), ...
    geodetic_frame_convention, use_WGS84);


%% Show the vehicle trajectory
plotTrajectory(traj, states_position, geodetic_frame_convention);


%% Generate and save sensor data
disp('Generating sensor data...')
for lv1 = 1 : num_samples
    % Generate data
    sensorData.meas_position = myPositionSensor.generateMeasurement(states_position.r_bt_t);
    [sensorData.meas_odom_lin, sensorData.meas_odom_ang] = ...
                        myOdomSensor.generateMeasurement(states_odom.v_btt_b, states_odom.omega_bt_b);
    
    % Include earth's rotation in IMU measurements
    if include_earth_rotation
        [sensorData.meas_accel, sensorData.meas_gyro, gt_accel_bias, gt_gyro_bias] = ...
                     myIMU.generateMeasurement(states_imu.f_b_in, states_imu.omega_bi_b);
    % Treat local geodetic frame as inertial frame
    else
        [sensorData.meas_accel, sensorData.meas_gyro, gt_accel_bias, gt_gyro_bias] = ...
                     myIMU.generateMeasurement(states_imu.f_b_tan, states_imu.omega_bt_b);
    end
    % Save data
    sensorData.(['gt_accel_bias_', num2str(lv1)]) = gt_accel_bias;
    sensorData.(['gt_gyro_bias_', num2str(lv1)]) = gt_gyro_bias;
    %save(['output/sensor_data_', num2str(lv1)], '-struct', 'sensorData');   
end


%% Save ground-truth states and trajectory info
% Save intero or extero states depending on preint
if IMU_preint
    gt_states = states_position;
else
    gt_states = states_imu;
end

save('output/gt_states', 'gt_states');
save('output/trajectory_config', 'traj')

% Assemble ground truth full state
gt_SE23R3R3 = [];
for idx = 1:length(gt_states.time)
    C_ab_i = gt_states.C_bt(:,:,idx)';
    v_a_i = gt_states.v_btt_t(:,idx);
    r_a_i = gt_states.r_bt_t(:,idx);
    b_g_i = sensorData.gt_gyro_bias_1(:,idx);
    b_a_i = sensorData.gt_accel_bias_1(:,idx);
    gt_SE23R3R3(:,:,idx) = SE23xR3x2.synthesize(C_ab_i, v_a_i, r_a_i, b_g_i, b_a_i);
end

%% Get sensor timestamps
IMU_f_inv = 1 / f_imu;
IMU_t_max = IMU_f_inv * (length(sensorData.meas_gyro) - 1);
sensorData.t_IMU = 0:IMU_f_inv:IMU_t_max;

pos_f_inv = 1 / f_position;
pos_t_max = pos_f_inv * (length(sensorData.meas_position) - 1);
sensorData.t_pos = 0:pos_f_inv:pos_t_max;

odom_f_inv = 1 / f_odom;
odom_t_max = odom_f_inv * (length(sensorData.meas_odom_lin) - 1);
sensorData.t_odom = 0:odom_f_inv:odom_t_max;

