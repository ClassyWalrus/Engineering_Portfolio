%========================================================================%
%                   X-15_NLPMEOM_Simulation.m
%
%   This m-script simulates the longitudinal point mass dynamics of the 
%   X-15 hypersonic aircraft by solving the nonlinear, longitudinal, point
%   mass, equations of motion (NLPMEOM). 
%
%   Inputs:
%
%       rho = density;
%       S = wing area
%       CD = drag coefficient at max L/D
%       CL = lift coefficient at max L/D
%       m = mass of aircraft (in kg)
%       x_0 = initial conditions (x = [U;G;d;h])
%
%   Outputs:
%
%       t_out = Time vector out
%       U = air speed magnitude
%       G = flight path angle
%       d = distance traveled
%       h = altitude
%
%   
%   Created:            February 2, 2024, Demoz Gebre-Egziabher
%   Last Modified:      February X, 2024 (I. M. Student)
%
%=========================================================================

%   (1)     Clear the workspace

close all;
clear variables;
clc;

%   (2)     Define Constants

%   (2.a)   Conversion and aircraft independent constants

ft2m = 0.3048;              %   Feet to meters
lb2kg = 0.4536;             %   Pounds to kilograms
r2d = 180/pi;               %   Radians to degrees
d2r = 1/r2d;                %   Radians to degrees
g = 9.81;                   %   Magnitude of Gravitational Acceleration

%   (2.b)   Atmopsheric Constants

rho = [0.653; 0.3025];      %   Density at FS1 and FS2
a = [316; 294.9];          %   Speed of sound at FS1 and FS2

%   (2.c)   X-15 aerodynamic and geometric constants

m = 15560*lb2kg;                 %   Aircraft mass in kg
S = 200*(ft2m)^2;                %   Wing planform area in meters square
CD = [0.686; 0.19806];                %   CD at max L/D for FS1 and FS2
CL = [0.332; 0.6519];                %   CL at max L/D for FS1 and FS2

%   (2.d)   Time vector and intial condition

dt = .1;
t_in = 0;
x_0 = [a(1)*0.5 0  0 20000*ft2m;...
        a(2)*1 0 0 40000*ft2m];


%   Define the EOM as an annonymous MATLAB function

FS = 2;                         %   Flight scenario flag

NLPMEOM = @(t,x) [-rho(FS)*x(1)^2*S*CD(FS)/(2*m) - g*sin(x(2));...
                  rho(FS)*x(1)*S*CL(FS)/(2*m) - g/x(1)*cos(x(2));...
                    x(1)*cos(x(2));...
                    x(1)*sin(x(2))];

%   Solve the NLPMEOM by calling ODE45

[t_out,x_out] = ode45(NLPMEOM, [0:dt:5*60] ,x_0(FS,:));

%   Package output into appropriately named variables for plotting

t = t_out;
U = x_out(:,1);
G = x_out(:,2)*r2d;
d = x_out(:,3)/1000;
h = x_out(:,4);

%   Plot results

figure
subplot(221)
plot(t, U,'r-','LineWidth',2);grid on;
ylabel('$U$ (m/s)','Interpreter','latex');
subplot(222)
plot(t,G,'r-','LineWidth',2);grid on;
ylabel('$\Gamma$ (deg)','Interpreter','latex');
subplot(223)
plot(t,d,'r-','LineWidth',2);grid on;
ylabel('$d$ (km)','Interpreter','latex');
xlabel('$t$ (sec)','Interpreter','latex');
subplot(224)
plot(t,h,'r-','LineWidth',2);grid on;
ylabel('$h$ (m)','Interpreter','latex');
xlabel('$t$ (sec)','Interpreter','latex');

