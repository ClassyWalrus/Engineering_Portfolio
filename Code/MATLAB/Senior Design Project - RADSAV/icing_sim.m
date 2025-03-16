% Icing Simulation with Pitot Tube (Rounded Nose) Geometry
clear; clc; close all;
% Vehicle and Atmospheric Parameters
v_inf = @(z) 10*exp(-z/1800); % Velocity as a function of altitude (m/s)
z_max = 40e3; % Maximum altitude (m)
t_total = 3600; % Total time (s)
dt = 1e-2; % Time step (s)
N_t = t_total / dt;
% Atmospheric model functions (simplified)
T_air_func = @(z) max(288.15 - 0.0065*z, 216.65); % Temperature (K)
P_air_func = @(z) 101325 * (1 - 0.0065*z/288.15)^5.256; % Pressure (Pa)
rho_air_func = @(z) P_air_func(z) / (287.05 * T_air_func(z)); % Density (kg/m^3)
LWC_func = @(z) max(0.1 - z/12000, 0); % Liquid water content (g/m^3)
h_conv_stag = @(z, v) 100 + 10*sqrt(rho_air_func(z)*v); % Stagnation heat transfer coefficient (W/m^2·K)
% Geometry
nose_radius = 0.01*(1/2.54); % Radius of the pitot tube nose (m)
body_length = 1; % Length of the cylindrical body (m)
dx = 0.01; % Spatial step (m)
% Discretize geometry
x_nose = 0:dx:(pi * nose_radius); % Nose region (semicircular)
x_body = (pi * nose_radius + dx):dx:(pi * nose_radius + body_length); % Body region
x = [x_nose, x_body]; % Combine regions
Nx = length(x);
% Initialize variables
ice_thickness = zeros(1, Nx); % Ice thickness (m)
surface_temp = 303.15 * ones(1, Nx); % Surface temperature (K)
altitude = linspace(0, z_max, N_t); % Altitude (m)
% Main loop
for t = 1:N_t
    z = altitude(t); % Current altitude (m)
    T_air = T_air_func(z);
    LWC = LWC_func(z);
    rho_air = rho_air_func(z);
    v = v_inf(z);
    for i = 1:Nx
        % Determine region (nose or body)
        if i <= length(x_nose)
            % Nose region (stagnation effects)
            theta = x_nose(i) / nose_radius; % Angular position
            h_conv = h_conv_stag(z, v) * cos(theta); % Heat transfer
            eta = 0.1 * cos(theta); % Collection efficiency decreases with angle
        else
            % Body region
            h_conv = 100 + 10*rho_air*v^0.8; % Simplified body heat transfer
            eta = 0.05; % Uniform collection efficiency
        end
        % Mass flux of droplets impacting the surface
        m_dot_droplet = LWC * v * eta;
        % Heat flux due to convection
        q_conv = h_conv * (T_air - surface_temp(i));
        % Heat flux required for freezing
        q_freeze = m_dot_droplet * 334e3; % Latent heat of fusion
        % Update ice thickness
        if q_conv < q_freeze
            ice_accumulation_rate = (q_freeze - q_conv) / (917 * 334e3);
            ice_thickness(i) = ice_thickness(i) + ice_accumulation_rate * dt;
        end
    end
    % Visualization
    if mod(t, 10) == 0
        plot(x, ice_thickness, 'b-', 'LineWidth', 2);
        title(sprintf('Ice Thickness at t = %.1f s, Altitude = %.1f km', t*dt, z/1000));
        xlabel('Position along the surface (m)');
        ylabel('Ice Thickness (m)');
        grid on;
        drawnow;
    end
end
% Final plot
figure;
plot(x, ice_thickness, 'r-', 'LineWidth', 2);
title('Final Ice Thickness Distribution');
xlabel('Position along the surface (m)');
ylabel('Ice Thickness (m)');
grid on;
