% AUTHOR:       Daniil Lisus (260669654)
% TOPIC:        MLG Batch Estimation
% MODULE:       Trajectory Generation
% DESCRIPTION:  

classdef robot_path

    properties (Access = private)
        
    end
    
    properties
        name;
        dt;             % [s]
        ideal_path;     % [T or F]
        theta;          % [rad]
        omega_true;     % [rad/s]
        omega_meas;     % [rad/s]
        sin_theta;
        cos_theta;
        r_true;         % [m, m, m]
        r_meas;         % [m, m, m]
        v_true;         % [m/s, m/s, m/s]
        v_meas;         % [m/s, m/s, m/s]
        t;              % [s]
        C_ab;           % DCM
        xi_phi;         % xi_phi from DCM
        b;              % [rad/s] bias
        q;              % quaternion [x y z w]
        r_uz_b;         % correction term for shifting position from pivot 
                        % point on the body (z) to the uwb center (u) in
                        % body frame (a)
        X_SE3;          % SE(3) state
        Xi_SE3;
    end
    
    methods
        function path = robot_path(name,struct_name,r_uz_b,ideal_path,t,r_0,theta_0,omega,v)
            path.name = name;
            if nargin >1
                path.r_uz_b = r_uz_b;

                if strcmp(struct_name,"None")
                    if ~isequal(size(omega),size(v))
                        error('Omega, v dimensions must be the same!')
                    end
                    path.t = t;
                    path.dt = t(2)-t(1);

                    path.ideal_path = ideal_path;
                    path.omega_true = omega;
                    path.v_true = v;

                    % Generate full path
                    path = path.generatePath(r_0, theta_0);

                    % Corrupt measurements if not ideal
                    if ~path.ideal_path
                        path = path.corruptMeasurements;
                    else
                        path.r_meas = path.r_true;
                        path.omega_meas = path.omega_true;
                        path.v_meas = v;
                        path.b = [0; 0; 0]*t;
                    end
                    for ii = 1:length(t)
                        path.Xi_SE3(:,:,ii) = SE3.synthesize(path.dt*path.omega_meas(:,ii), path.dt*path.v_meas(:,ii));
                    end
                else
                    path.t = struct_name.t;

                    for ii = 1:length(struct_name.t)
                        C_ab = struct_name.C_ba(:,:,ii)';
                        C_ab = (C_ab*C_ab')^(-0.5)*C_ab;
                        phi = SO3.decompose(C_ab);
                        path.xi_phi(:,ii) = phi;
                        path.theta(ii,:) = phi(3);
                        path.sin_theta(ii,:) = C_ab(2,1);
                        path.cos_theta(ii,:) = C_ab(1,1);
                        path.C_ab(:,:,ii) = C_ab;
                        path.r_meas(:,ii) = struct_name.r_zw_a(:,ii) + C_ab*r_uz_b;
                    end

                    q_temp = dcmToQuat(path.C_ab)';
                    path.q = [q_temp(:,2) q_temp(:,3) q_temp(:,4) q_temp(:,1)];

                    % Get rid of nan values
                    for ii = 2:length(path.r_meas(1,:))
                       if any(isnan(path.r_meas(:,ii)))
                           path.r_meas(:,ii) = path.r_meas(:,ii-1);
                       end
                    end
                    for ii = 2:length(path.theta)
                       if isnan(path.theta(ii))
                           path.theta(ii) = path.theta(ii-1);
                       end
                    end
                    path.r_true = path.r_meas;
                    path.omega_meas = struct_name.gyro;
                    path.b = 0.0*path.r_meas;
                end
            end
        end
        
        function ori = getOrientation(path, timestamp, zero_start)
            time_step = int32(timestamp/path.dt) + zero_start;     % Transform timestamp to searchable index
            ori = path.theta(time_step);
        end
        
        function pos = getPosition(path, timestamp, zero_start)
            time_step = int32(timestamp/path.dt) + zero_start;     % Transform timestamp to searchable index
            pos = [path.r_meas(1,time_step); path.r_meas(2,time_step); path.r_meas(3,time_step)];
        end
        
        function dist = getRange(this_path,other_path, timestamp, zero_start)
            parameters;
            pos1 = this_path.getPosition(timestamp, zero_start);
            pos2 = other_path.getPosition(timestamp, zero_start);
            
            dist = norm(pos1-pos2);
            
            if ~this_path.ideal_path
                % Generate noise
                r_noise = sqrt(range_var).*randn(1)';     % Measurement x pos noise
                dist = dist + r_noise;
            end
        end

    end
    
    methods (Access = private)
        function path = generatePath(path, r_0, theta_0)
            theta_ideal(1) = theta_0(3);
            sin_ideal(1) = sin(theta_ideal(1));
            cos_ideal(1) = cos(theta_ideal(1));
            r_ideal(:,1) = r_0;
            
            C_0 = SO3.synthesize(theta_0);
            T_0 = SE3.synthesize(C_0, r_0);
            C_ideal(:,:,1) = C_0;
            xi_phi(:,1) = theta_0;
            
            T_k = T_0;
            for kk = 2:length(path.t)
                T_k_min = T_k;

                % Get previous velocities
                u1_k_min = path.omega_true(:, kk-1);
                u2_k_min = path.v_true(:, kk-1);

                % Calculate full Xi matrix
                psi = so3alg.expMap(so3alg.wedge(path.dt * u1_k_min));
                dtu2 = path.dt*u2_k_min;
                Xi_k_min = [psi dtu2; zeros(1,3) 1];

                % Get updated T matrix
                T_k = T_k_min*Xi_k_min;

                % Store values
                [ C_k, r_k ] = SE3.decompose(T_k);
                angles = SO3.decompose(C_k);
                theta_ideal(kk) = angles(3);
                sin_ideal(kk) = C_k(2,1);
                cos_ideal(kk) = C_k(1,1);
                r_ideal(:, kk) = r_k;
                C_ideal(:,:,kk) = C_k;
                xi_phi(:, kk) = angles;
            end
            
            % Save path
            path.theta = theta_ideal;
            path.r_true = r_ideal;
            path.sin_theta = sin_ideal; 
            path.cos_theta = cos_ideal;
            path.C_ab = C_ideal;
            path.xi_phi = xi_phi;
            %path.X_SE3(:,:,1:size(r_ideal,2)) = eye(4);
            path.X_SE3(1:3,1:3,:) = C_ideal;
            path.X_SE3(1:3,4,:)= r_ideal;
            path.X_SE3(4,4,:)= ones(1,1,size(r_ideal,2));
        end
        
        function path = corruptMeasurements(path) 
            % Generate noise
            parameters;
            % GPS
            v_r1 = GPS_x_std.*randn(size(path.t));    % Measurement x pos noise
            v_r2 = GPS_y_std.*randn(size(path.t));    % Measurement y pos noise
            v_r3 = GPS_z_std.*randn(size(path.t));    % Measurement y pos noise
            % Gyro
            w_o1 = gyro_x_std.*randn(size(path.t));         % Measurement x gyro noise
            w_o2 = gyro_y_std.*randn(size(path.t));         % Measurement y gyro noise
            w_o3 = gyro_z_std.*randn(size(path.t));         % Measurement y gyro noise
            g_b1 = gyro_x_bias + path.dt*gyro_x_walk.*randn(size(path.t));
            g_b2 = gyro_y_bias + path.dt*gyro_y_walk.*randn(size(path.t));
            g_b3 = gyro_z_bias + path.dt*gyro_z_walk.*randn(size(path.t));
            % Velocity
            w_v1 = vel_x_std.*randn(size(path.t));          % Measurement x vel noise
            w_v2 = vel_y_std.*randn(size(path.t));          % Measurement y vel noise
            w_v3 = vel_z_std.*randn(size(path.t));          % Measurement z vel noise
            
            % Populate noise
            r_cor = path.r_true - [v_r1; v_r2; v_r3];
            omega_cor = path.omega_true - [w_o1; w_o2; w_o3] - [g_b1; g_b2; g_b3];
            vel_cor = path.v_true - [w_v1; w_v2; w_v3];
            
            % Save path
            path.r_meas = r_cor;
            path.omega_meas = omega_cor;
            path.b = [g_b1; g_b2; g_b3];
            path.v_meas = vel_cor;
        end
        
    end
    
    methods (Static)
        
        function plotPaths(paths)
            parameters;     % Need to port over plotting parameters
            
            figure
            hold on
            grid on
            title('\textbf{Robot Path(s)}','fontsize',font_size,'Interpreter','latex') 
            xlabel('$x$ (m)','fontsize',font_size,'Interpreter','latex');
            ylabel('$y$ (m)','fontsize',font_size,'Interpreter','latex');
            max_x = 0;
            min_x = 0;
            max_y = 0;
            min_y = 0;
            ii = 1;
            for path = paths
                start_pt = [path.r_true(1,1), path.r_true(2,1)];
                end_pt = [path.r_true(1,end), path.r_true(2,end)];
                h0(ii) = plot(path.r_true(1,:),path.r_true(2,:),'o','Linewidth',line_width,'DisplayName',sprintf('%s', path.name));
                %h1 = plot(start_pt(1),start_pt(2),'g*','Linewidth',line_width,'DisplayName','Start Points');
                %h2 = plot(end_pt(1),end_pt(2),'r*','Linewidth',line_width,'DisplayName','End Points');
                txt = sprintf('%s \\rightarrow ', path.name);
                text(start_pt(1),start_pt(2),txt,'HorizontalAlignment','right')
                max_x = max(max(path.r_meas(1,:)),max_x);
                min_x = min(min(path.r_meas(1,:)),min_x);
                max_y = max(max(path.r_meas(2,:)),max_y);
                min_y = min(min(path.r_meas(2,:)),min_y);
                ii = ii+1;
            end
            xlim([min_x-5 max_x+5])
            ylim([min_y-5 max_y+5])
            legend([h0],'Interpreter', 'latex')
            set(gca,'FontSize',font_size_legend);
        end
        
        function plotOrientations(paths)
            parameters;     % Need to port over plotting parameters
            
            figure
            hold on
            grid on
            title('\textbf{Robot Orientation(s)}','fontsize',font_size,'Interpreter','latex') 
            xlabel('$t$ (m)','fontsize',font_size,'Interpreter','latex');
            ylabel('$\theta$ (m)','fontsize',font_size,'Interpreter','latex');
            max_t = 0;
            min_t = 0;
            max_theta = 0;
            min_theta = 0;
            ii = 1;
            for path = paths
                t = path.t;
                theta = path.theta;
                start_pt = [t(1), theta(1)];
                end_pt = [t(end), theta(end)];
                h0(ii) = plot(t,theta,'o','Linewidth',line_width,'DisplayName',sprintf('%s', path.name));
                h1 = plot(start_pt(1),start_pt(2),'g*','Linewidth',line_width,'DisplayName','Start Points');
                h2 = plot(end_pt(1),end_pt(2),'r*','Linewidth',line_width,'DisplayName','End Points');
                max_t = max(max(t),max_t);
                min_t = min(min(t),min_t);
                max_theta = max(max(theta),max_theta);
                min_theta = min(min(theta),min_theta);
                ii = ii+1;
            end
            xlim([min_t-5 max_t+5])
            ylim([min_theta-5 max_theta+5])
            legend([h0, h1(1), h2(1)],'Interpreter', 'latex')
            set(gca,'FontSize',font_size_legend);
        end
   end
end

