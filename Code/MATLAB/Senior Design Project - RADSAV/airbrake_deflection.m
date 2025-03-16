%% Avery Herdina
%% Senior Design
%% Airbrake Deformation Deflection

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Assumed standard atmosphere conditions:
% Altitude: 30,000m
% Mach: 1.5

close all;
clear all;

% Convert altitude to meters and calculate atmospheric properties
altitude = 30000; %altitude in meters
[T, a, p, rho] = atmosisa(altitude, 'extended', 'on'); % Obtain standard atmospheric parameters

% Mach and Velocity
M = 1.5; % Mach Number
y = 1.4; %gamma
V = a * M; % Velocity (m/s)
P_stagnation = ((((1+((1.4-1)/2)*M^2))^((-1.4)/(1.4-1)))/p)^-1; %(Pa)
T_stagnation = (((1+((1.4-1)/2)*M^2)^-1)/T)^-1;
%q = 1.4/2*p*M^2; %Dynamic pressure in isentropic flow used for force on airbrake
%q = p*(1+((1.4-1)/2)*M^2)^(1.4/(1.4-1))-p;

q = p*( (((((y + 1)/2)*(M^2))^(y/(y-1))) / (((2*y*M^2)/(y + 1)) - (y-1)/(y+1))^(1/(y-1))) - 1);

% Force on airbrake
Cd = 2.21; % Coefficient of drag
A = 0.0029; % Surface area of 1 exposed airbrake leaf in m^2
F = Cd * A * q; % Force on Airbrake (N)
P = F / A; % Pressure on Airbrake (Pa)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Cross-section and moment of inertia for beam deflection
b = convlength(4.3, 'in', 'm'); % Width of cross-section (m)
t_in = [0.5, 0.25, 0.125];
t = convlength([0.5, 0.25, 0.125], 'in', 'm'); % Thicknesses of airbrake (m)
I = (b .* t.^3) ./ 12; % Moment of inertia (m^4)

% Deflection of Airbrake treated as a cantilever beam
L = linspace(0, convlength(1.07, 'in', 'm'), 1000); % Fully extended length of airbrake
E = 29E9; % Elastic modulus (Pa)

% Calculate deflection for each thickness
delta = zeros(length(t), length(L));
for i = 1:length(t)
    delta(i, :) = (F .* L.^3) ./ (3.* E .* I(i)); % Deflection formula
end

% Find the maximum deflection across all thicknesses
max_deflection = max(delta, [], 2);
global_y_max = max(max_deflection); % Global max deflection for y-axis scaling

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Plotting deflections
figure;

for i = 1:length(t)
    subplot(3, 1, i);
    plot(L, delta(i, :), 'LineWidth', 1.5);
    hold on;
    
    % Mark and label the maximum deflection
    [max_val, max_idx] = max(delta(i, :));
    plot(L(max_idx), max_val, 'ro', 'MarkerSize', 8, 'LineWidth', 1.5); % Red marker
    text(L(max_idx), max_val, sprintf('Max: %.2e m', max_val), ...
        'VerticalAlignment', 'bottom', 'HorizontalAlignment', 'right', ...
        'FontSize', 10, 'Color', 'r');
    
    % Title and labels
    title(sprintf('Deflection of Airbrake (%.3f in Thick)', t_in(i)), 'FontSize', 12);
    xlabel('Length (m)');
    ylabel('Deflection (m)');
    grid on;

    % Set the y-axis to the same range for all subplots
    ylim([0, 1.1 * global_y_max]);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
disp('Calculation Complete!');
