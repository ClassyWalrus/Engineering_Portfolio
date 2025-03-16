close all;
clear all
clc

%% Mechanics of Solids HW5
% Colton Davies

%% Question 1:

%Initial Params
Nx1 = 100;
Nx2 = 0.2*Nx1;
k = 0.02;
x1 = 0:1:Nx1;
x2 = -Nx2:1:Nx2;
[X1,X2] = meshgrid(x1,x2);

%Plot Undeformed Configuration
figure(1)
subplot(1,2,1)
scatter(X1(:),X2(:), 4,'b');
grid on;
title('Reference Configuration');
xlabel('$\tilde{e}_{1}$', 'Interpreter', 'latex');
ylabel('$\tilde{e}_{2}$', 'Interpreter', 'latex');
axis equal

%Deformations
y1 = (1/k)*sin(k*X1) - X2.*sin(k*X1);
y2 = (-1/k)*(cos(k*X1)-1) + X2.*cos(k*X1);

%Deformation Gradient
F = [cos(k*X1) - k*X2.*cos(k*X1), -sin(k*X1);
     sin(k*X1) - k*X2.*sin(k*X1), -cos(k*X2)];


%Plot Deformed Configuration
subplot(1,2,2)
scatter(y1(:), y2(:), 4, 'r')
grid on;
title('Deformed Configuration');
xlabel('$\tilde{e}_{1}$', 'Interpreter', 'latex');
ylabel('$\tilde{e}_{2}$', 'Interpreter', 'latex');
axis equal

%Attempted to correlate the magnitude of deformation to a color gradient,
%but I was dedicating too much time to learning the syntax vs learning the
%content


%% Question 2:

%Initial Params
k = 0.02;
Epsilon_Vals = [1.5, 1, 0.5];
i = -5:1:5;
j = -5:1:5;
X_0 = [75, 0]; % Point of Interest in the Body

figure(2)
for g = 1:length(Epsilon_Vals)
    
    %Find Square Area in Body
    E = Epsilon_Vals(g);
    [I, J] = meshgrid(i, j);
    X1_Small = X_0(1) + E * I;
    X2_Small = X_0(2) + E * J;

    %Deformations
    y1 = (1/k)*sin(k*X1_Small) - X2_Small.*sin(k*X1_Small);
    y2 = (-1/k)*(cos(k*X1_Small)-1) + X2_Small.*cos(k*X1_Small);

    % Plot Deformations for each Epsilon
    subplot(1, length(Epsilon_Vals), g);
    scatter(y1(:), y2(:), 4, 'r');
    title(['\epsilon = ', num2str(E)]);
    xlabel('$\tilde{e}_{1}$', 'Interpreter', 'latex');
    ylabel('$\tilde{e}_{2}$', 'Interpreter', 'latex');
    xlim([40, 60]);
    ylim([35, 55]);
    grid on
end

%% Question 3:

% Initial Constants
Kappa = 0.6; % Shear 
Delta = 0.2; % Extension 
Epsilon = 0.2; % Extension 

% Simple Shear Deformation
y1_Shear = X1 + Kappa * X2;
y2_Shear = X2;

% Simple Extension Deformation
y1_Extension = (1 - Delta) * X1;
y2_Extension = (1 + Epsilon) * X2;

% Plot Undeformed, Simple Shear, and Simple Extension
figure(3)
subplot(1,3,1)
scatter(X1(:), X2(:), 6, 'b');
title('Undeformed Configuration');
xlabel('$\tilde{e}_{1}$', 'Interpreter', 'latex');
ylabel('$\tilde{e}_{1}$', 'Interpreter', 'latex');
xlim([-20, 120]);
ylim([-50, 50]);
grid on

subplot(1,3,2)
scatter(y1_Shear(:), y2_Shear(:), 6, 'r');
title('Simple Shear');
xlabel('$\tilde{e}_{1}$', 'Interpreter', 'latex');
ylabel('$\tilde{e}_{1}$', 'Interpreter', 'latex');
xlim([-20, 120]);
ylim([-50, 50]);
grid on

subplot(1,3,3)
scatter(y1_Extension(:), y2_Extension(:), 6, 'k');
title('Simple Extension');
xlabel('$\tilde{e}_{1}$', 'Interpreter', 'latex');
ylabel('$\tilde{e}_{1}$', 'Interpreter', 'latex');
xlim([-20, 120]);
ylim([-50, 50]);
grid on


%% Question 4:

% Initial Constants
Kappa = 0.8; % Shear constant

% Deformation Gradient for Simple Shear
F = [1, Kappa; 0, 1]

% Left Cauchy-Green Strain Matrix B
B = F * F'

% Eigenvalues and eigenvectors of B
[Eig_Vec, Eig_Val] = eig(B);
Eig_Val
e_1 = Eig_Vec(:,1)
e_2 = Eig_Vec(:,2)

% Plot a circle of points (radius = 1) centered at the origin
theta = linspace(0, 2*pi, 100);
X1_Circ = cos(theta);
Y1_Circ = sin(theta);

Circ_Deform = F * [X1_Circ; Y1_Circ];
Eig_Deform = F * Eig_Vec;

% Plot the original circle
figure(4)
hold on;
scatter(X1_Circ(:), Y1_Circ(:), 10, 'b', 'filled');
scatter(Eig_Vec(1,1), Eig_Vec(2,1), 35, 'g', 'filled');
scatter(Eig_Vec(1,2), Eig_Vec(2,2), 35, 'k', 'filled');
axis equal;

% Plot the deformed circle
scatter(Circ_Deform(1,:), Circ_Deform(2,:), 10, 'r', 'filled');
scatter(Eig_Deform(1,1), Eig_Deform(2,1), 35, 'g', 'filled');
scatter(Eig_Deform(1,2), Eig_Deform(2,2), 35, 'k', 'filled');
title('Simple Shear of a Circle');
xlabel('$\tilde{e}_{1}$', 'Interpreter', 'latex');
ylabel('$\tilde{e}_{1}$', 'Interpreter', 'latex')
grid on