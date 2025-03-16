clear all
close all

%% Mechanics of Solids HW4

%Define material params
L = 1; %m
E = 2e9; %Pa
I = (pi/4)*(0.01^4); %m^4
R_cr = (pi^2*E*I) / (4*L^2); %N
ns = 1000;
x_0 = 0;
y_0 = 0;
N= 6;  % Number of beams and division of R_cr

% Pre-allocate deflections
Deflection_x = zeros(N, ns);
Deflection_y = zeros(N, ns);

for i = 1:N
    %Compute load for each i
    R = (1+(i/N))*R_cr;
    c = sqrt((2*R) / (E*I));

    %Use elliptic integral to find phi(x)
    Phi_s = @(phi) ellipticK(sin(phi)^2) - (c*L/sqrt(2));
    Phi_L = fzero(Phi_s, pi/4);

    %Use elliptic integral to find s(X) and invert for X(s)
    X_vec = linspace(0, pi/2, ns);
    s_X = ellipticF(X_vec, sin(Phi_L)^2) * (sqrt(2)/c);
    % [Phi_0,Phi_L] = ellipke(M); need to break into complete + incomplete
    X_s = [s_X; X_vec];

    % Compute the angle theta(s) (related to the end of the beam)
    % theta_s = 2*Phi_s;
    theta_s = 2*asin(sin(Phi_L) * sin(X_s(2,:)));

    %Take integrals
    int_x = cumtrapz(s_X, cos(theta_s));
    int_y = cumtrapz(s_X, sin(theta_s));

    %Deflections
    U_x = x_0 + int_x;
    U_y = y_0 + int_y;
Deflection_x(i,:) = U_x;
Deflection_y(i,:) = abs(U_y);
end

%Graphing
figure
hold on
xlabel('y_1');
ylabel('y_2');
title('Beam Deflection', 'FontSize', 14);
for j = 1:i
    plot(Deflection_x(j,:), Deflection_y(j,:), 'LineWidth', 2)
end
grid on;
legend(arrayfun(@(x) sprintf('R = (1+1/%.1f) * R_{cr}', x), (1:N)', 'UniformOutput', false));
legend show;
hold off;