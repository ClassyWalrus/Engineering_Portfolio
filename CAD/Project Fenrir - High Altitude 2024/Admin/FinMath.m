% Code for determining the Leading and Trailing edge chamfer 
% length
%
% Daniel Cockcroft
% Rocket Team

% INPUTS
t = 0.21;                           % Thickness (inches)
r = 0.02;                           % Radius of leading/trailing edge of edge of fin (inches)
Chamfer_Angle = 15;                 % Angle (degrees)
X = (t/2-r)/tand(Chamfer_Angle);    % Calculate the length the bevel cuts at Gary's (inches)

%% SAC 2024 Geometry
% L_Sweep = 10;       % Sweep length
% L_Tip = 5;          % Tip chord length
% H = 6.75;           % Fin height
% L_Root = 16;

%% HA 2024 Fin Geometry
% Sustainer Fins
% L_Sweep = 10;       % Sweep length
% L_Tip = 3;          % Tip chord length
% H = 5.5;            % Fin height
% L_Root = 12;        % Root chord

%% Booster Fins Geometry
L_Sweep = 11;         % Sweep length
L_Tip = 4;            % Tip chord length
H = 6;                % Fin height
L_Root = 16;          % Root chord


% Angles of leading and trailing edge from horizontal
theta = atand(H/L_Sweep);
alpha = abs(atand(H/(L_Root-(L_Sweep+L_Tip))));


% Leading and Trailing edge input into RASAero
X_LE = X/cosd(90-theta)
X_TE = X/cosd(90-alpha)