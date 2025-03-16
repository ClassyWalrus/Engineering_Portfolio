% Created by Prof. Caverly
% Updated Fall 2024

% This is a sample code for AEM 4602W that is used to solve the 
% 2-dimensional GPS navigation solution without any clock errors

clear all       % clears the workspace
clc             % clears the command window
close all       % closes all figures

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Define Constants

% Position of the pseudolites (TO BE ENTERED)
r_S1 = [2.03;3.79];
r_S2 = [-17.56;9.25];
r_S3 = [4.40;8.61];

% Ranges between the pseudolites and the receiver (TO BE ENTERED)
rho1 = 3.49;
rho2 = 16.95;
rho3 = 6.29;

% Initial guess of receiver position
r_R_bar = [0;0];

% Maximum number of iterations for Newton's method
max_iters = 10;

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Solve Using Newton's Method

% Initialized loop variables
iters = 0;
stop = 0;
while stop == 0
	iters = iters+1; % update number of iterations

	% Compute range based on current estimate of receiver position
	rho1_bar = norm(r_R_bar-r_S1);
	rho2_bar = norm(r_R_bar-r_S2);
	rho3_bar = norm(r_R_bar-r_S3);

	% Compute matrices from linearization for Newton's method
	A = [(r_R_bar(1)-r_S1(1))/rho1_bar (r_R_bar(2)-r_S1(2))/rho1_bar;...
		(r_R_bar(1)-r_S2(1))/rho2_bar (r_R_bar(2)-r_S2(2))/rho2_bar;...
		(r_R_bar(1)-r_S3(1))/rho3_bar (r_R_bar(2)-r_S3(2))/rho3_bar];
	b = [rho1 - rho1_bar; rho2 - rho2_bar; rho3 - rho3_bar];

	% Solve for updated receiver position estimate
	delta_r = A\b;
	r_R_bar = r_R_bar + delta_r;

	% Stopping criteria based on convergence on solution within 1e-6 m or
	% exceeding maximum number of iterations
	if ((norm(delta_r) < 1e-6) || (iters > max_iters)) 
		stop = 1;
	end
end

% Estimated position of the receiver
r_R_sol = r_R_bar



