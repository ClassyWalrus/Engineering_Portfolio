% Created by Prof. Caverly
% Updated Fall 2024
% This is a sample code for AEM 4602W that is used to solve the 
% 2-dimensional GPS navigation solution without any clock errors
clear all       % clears the workspace
clc             % clears the command window
close all       % closes all figures
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Import Values %%

cal = ["cal_data0_Extra2.mat", "cal_data1_Extra2.mat", "cal_data2_Extra2.mat", "cal_data3_Extra2.mat"];

for j = 1:4

load(cal(j));

%% Define Constants
% Extract satellite positions and pseudorange data
range_data = pseudo_range;           % Nx6 matrix for range measurements across time steps
average_pos = mean(range_data);
for i = 1:6
    err(i) = std(range_data(:,i));
end
PseudoStd = mean(err);

% individual_errors(j) = sum(err)/length(err);
end
% error_sum = sum(individual_errors)/length(individual_errors);

Flight = ["flight_5B_M_1_Extra2.mat", "flight_5B_M_2_Extra2.mat"];
for k = 1:2
%% Flight Data
load(Flight(k));


%% Define Constants
% Extract satellite positions and pseudorange data
sat_positions = GPS_sat_loc;        % 3x6 matrix for satellite positions
range_data = pseudo_range;           % Nx6 matrix for range measurements across time steps
num_timesteps = size(range_data, 1); % Total number of time steps
% Initialize variables to store estimated 3D positions for each time step
drone_positions = zeros(num_timesteps, 3);
drone_error = zeros(num_timesteps, 3);
% Initial guess for receiver position in 3D
r_R_bar = [0; 0; 0];   % Initial guess
max_iters = 100;        % Maximum iterations for Newton's method
% Loop over each time step to calculate 3D position
for t = 1:num_timesteps
    stop = 0;
    iters = 0;
    range_t = range_data(t, :)';  % Range data for current time step
    while stop == 0
        iters = iters + 1;
        
        % Compute estimated ranges for current position estimate
        rho_bar = zeros(6, 1);
        for i = 1:6
            rho_bar(i) = norm(r_R_bar - sat_positions(:, i)); % Range to each satellite
        end
        
        % Residuals (difference between actual and estimated ranges)
        delta_rho = range_t - rho_bar;
        
        % Construct the Jacobian matrix for 6 pseudolites in 3D
        J = zeros(6, 4);
        for i = 1:6
             % J(i, :) = [(r_R_bar' - sat_positions(:, i)') / rho_bar(i) 1];
             J(i, 1:3) = (r_R_bar' - sat_positions(:, i)') / rho_bar(i);
             J(i, 4) = 1;
        end
        
        % Update the receiver position estimate
        delta_r_R = J(:, 1:3) \ delta_rho;
        r_R_bar = r_R_bar + delta_r_R;
        
        % Check for convergence
        if norm(delta_r_R) < 1e-6 || iters >= max_iters
            stop = 1;
        end
    end
    
    % Store the estimated position for this time step
    drone_positions(t, :) = r_R_bar';
% Position Error calcs
q_mat = inv(J'*J);
p_dop = sqrt(q_mat(1,1)+q_mat(2,2)+q_mat(3,3));
t_dop = sqrt(q_mat(4,4));
n_dop = sqrt(q_mat(1,1));
e_dop = sqrt(q_mat(2,2));
d_dop = sqrt(q_mat(3,3));
drone_error(t, 1) = n_dop*3*PseudoStd;
drone_error(t, 2) = e_dop*3*PseudoStd;
drone_error(t, 3) = d_dop*3*PseudoStd;


end
%%
% Plotting the 3D trajectory of the drone
figure;
plot3(drone_positions(:, 1), drone_positions(:, 2) - 2.5, drone_positions(:, 3)- 3.8, 'k-', 'LineWidth', 1.5);
hold on;
% Plot satellite positions as blue circles
plot3(sat_positions(1, :), sat_positions(2, :), sat_positions(3, :), 'bo', 'MarkerSize', 8, 'DisplayName', 'Pseudolites');
% Error Bars
plot3(drone_positions(:,1) + drone_error(:,1), ...
      drone_positions(:,2) + drone_error(:,2) - 2.5, ...
      drone_positions(:,3) + drone_error(:,3) - 3.8, 'r--', 'LineWidth', 1.5);
plot3(drone_positions(:,1) - drone_error(:,1), ...
      drone_positions(:,2) - drone_error(:,2) - 2.5, ...
      drone_positions(:,3) - drone_error(:,3) - 3.8, 'r--', 'LineWidth', 1.5);
% Customize plot
xlabel('X Position (m)');
ylabel('Y Position (m)');
zlabel('Z Position (m)');
title('3D Trajectory of Drone');
legend('Drone Path', 'Pseudolites','Flight ');
grid on;
axis equal;
view(3);
hold off;

figure
% Error Plot
t_plot = linspace(0,1,num_timesteps);
subplot(3,1,1)
plot(time_data, drone_error(:,1),'k')
grid on
ylabel("North Axis Error (m)")
xlabel("Flight Time (s)")
% title("Drone Position Error (95% Confidence Interval)")
fontsize(15, "points")
subplot(3,1,2)
plot(time_data, drone_error(:,2),'k')
grid on
ylabel("East Axis Error (m)")
xlabel("Flight Time (s)")
fontsize(15, "points")
subplot(3,1,3)
plot(time_data, drone_error(:,3),'k')
grid on
ylabel("Down Axis Error (m)")
xlabel("Flight Time (s)")
fontsize(15, "points")
end