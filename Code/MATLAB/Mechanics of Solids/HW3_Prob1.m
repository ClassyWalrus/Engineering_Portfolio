%% Mehanics of Solids HW3
%Colton Davies
clear all;
close all;

%% Problem 1
% Constants chosen based off of 6061 Aluminum
L = 1; %m
L_10 = 10; %m
theta_L = 1000; %K
theta_0 = 278; %K
gamma = 20; %W/m^2*K
kappa = 237; %W/m*K
rho_0 = 2.7; %kg/m^3

zeta = sqrt(rho_0*gamma/kappa);
x = linspace(0,L,100);
x_10 = linspace(0,L_10,100);

%Solve temp for each bar length
theta_X_1 = (theta_L-theta_0)/(exp(zeta*L)-exp(-zeta*L)) * (exp(zeta*x)-exp(-zeta*x)) + theta_0;
theta_X_10 = (theta_L-theta_0)/(exp(zeta*L_10)-exp(-zeta*L_10)) * (exp(zeta*x_10)-exp(-zeta*x_10)) + theta_0;

%Graphing
figure(1)
plot(x,theta_X_1,'k-', 'LineWidth', 2)
grid on
title('1m Long Bar')
ylabel('Temperature (K)')
xlabel('Position (m)')

figure(2)
plot(x_10,theta_X_10, 'LineWidth', 2)
grid on
title('10m Long Bar')
ylabel('Temperature (K)')
xlabel('Position (m)')

%% Problem 3
%Define conditions
Curve_L = 10;
x_0 = 0;
y_0 = 0;
x_02 = -2;
y_02 = -3;

%Define thetas
ns = 1000;
s = linspace(0, Curve_L, ns);
theta = pi*sin(s);
theta_02 = pi/2*sin(s).^2;

%Take integrals
int_x = cumtrapz(s, cos(theta));
int_y = cumtrapz(s, sin(theta));
int_x_02 = cumtrapz(s, cos(theta_02));
int_y_02 = cumtrapz(s, sin(theta_02));

x = x_0 + int_x;
y = y_0 + int_y;

x_2 = x_02 + int_x_02;
y_2 = y_02 + int_y_02;

%Graphing
figure(3)
hold on
plot(x, y, 'k-', 'LineWidth', 2)
plot(x_2, y_2, 'r-', 'LineWidth', 2)
xlabel('X');
ylabel('Y');
title('Plane Inextensible Curves', 'FontSize', 14);
grid on;
legend('U(s) x_0=0 y_0=0', 'U(s) x_0=-2 y_0=-3');
hold off
