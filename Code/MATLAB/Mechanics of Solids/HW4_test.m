clear all
close all

%% Material Properties
L = 1;            % Beam length (m)
E = 2e9;          % Young's Modulus (Pa)
I = (pi/4) * (0.01^4); % Moment of inertia (m^4)
R_cr = (pi^2 * E * I) / (4 * L^2); % Critical load (N)
ns = 1000;         % Number of points for integration
N = 6;             % Number of load cases

% Pre-allocate deflections
Deflection_x = zeros(N, ns);
Deflection_y = zeros(N, ns);

%% Main Loop over N Load Cases
for i = 1:N
    % Compute load for current case
    R = (1 + (i / N)) * R_cr;
    c = sqrt((2 * R) / (E * I));
    
    % Define the function for solving Phi_L
    Phi_s = @(phi) ellipticK(sin(phi)^2) - (c * L / sqrt(2));
    Phi_L = fzero(Phi_s, pi / 4);
    
    % Compute s(X) and invert for X(s)
    X_vec = linspace(0, pi / 2, ns);
    s_X = ellipticF(X_vec, sin(Phi_L)^2) * (sqrt(2) / c);
    
    % Compute the angle theta(s) (related to the end of the beam)
    theta_s = 2 * asin(sin(Phi_L) * sin(X_vec));
    
    % Integrate to find deflections
    int_x = cumtrapz(s_X, cos(theta_s)); % x-deflection
    int_y = cumtrapz(s_X, sin(theta_s)); % y-deflection

    % Store the deflections for plotting
    Deflection_x(i, :) = int_x;
    Deflection_y(i, :) = abs(int_y);
end

%% Plotting
figure;
hold on;
xlabel('y_1');
ylabel('y_2');
title('Beam Deflection', 'FontSize', 14);
% Plot each beam curve
for j = 1:N
    plot(Deflection_x(j, :), Deflection_y(j, :), 'LineWidth', 2);
end
grid on;
legend(arrayfun(@(x) sprintf('R = %.1f * R_{cr}', x), (1:N)', 'UniformOutput', false));
legend show;
hold off;