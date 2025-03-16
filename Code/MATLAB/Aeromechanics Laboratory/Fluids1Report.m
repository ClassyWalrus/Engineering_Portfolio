 % This was a record of some of the live calculations and checks done on the
% FLUIDS 1 Lab

%Data pulling code created by Simeon Shaffar
%Everything else created + edited by Colton Davies

clc; clear; close all;

%% Immediate Measurements

% NACA 0016 (symmetrical airfoil)
span_imperial = 27.25; % inches
chord_imperial  = 11; % inches

span = span_imperial / 39.37; % m
chord = chord_imperial / 39.37; % m

A = span*chord;

mass = 1.807; %kg

% Constants
rho = 1.226; % kg/m^3


%% Loading Actual Data
% Load the data from each set (and each AoA)


%20 m/s velocity and AoA sweep from -12 to 24 degrees
Vel20_AoAneg12 = readtable("Vel20_AoAneg12.csv");
Vel20_AoAneg9  = readtable("Vel20_AoAneg9.csv");
Vel20_AoAneg6  = readtable("Vel20_AoAneg6.csv");
Vel20_AoAneg3  = readtable("Vel20_AoAneg3.csv");
Vel20_AoA0     = readtable("Vel20_AoA0.csv");
Vel20_AoApos3  = readtable("Vel20_AoApos3.csv");
Vel20_AoApos6  = readtable("Vel20_AoApos6.csv");
Vel20_AoApos9  = readtable("Vel20_AoApos9.csv");
Vel20_AoApos12 = readtable("Vel20_AoApos12.csv");
Vel20_AoApos15 = readtable("Vel20_AoApos15.csv");
Vel20_AoApos18 = readtable("Vel20_AoApos18.csv");
Vel20_AoApos21 = readtable("Vel20_AoApos21.csv");
Vel20_AoApos24 = readtable("Vel20_AoApos24.csv");
%30 m/s velocity and AoA sweep from 0 to 24 degrees
Vel30_AoApos0  = readtable("Vel30_AoApos0.csv");
Vel30_AoApos3  = readtable("Vel30_AoApos3.csv");
Vel30_AoApos6  = readtable("Vel30_AoApos6.csv");
Vel30_AoApos9  = readtable("Vel30_AoApos9.csv");
Vel30_AoApos12 = readtable("Vel30_AoApos12.csv");
Vel30_AoApos15 = readtable("Vel30_AoApos15.csv");
Vel30_AoApos18 = readtable("Vel30_AoApos18.csv");
Vel30_AoApos21 = readtable("Vel30_AoApos21.csv");
Vel30_AoApos24 = readtable("Vel30_AoApos24.csv");


% Calculate the mean of each dataset
avg_Vel20_AoAneg12 = mean( Vel20_AoAneg12{2:end, 2:end-1} , 1);
avg_Vel20_AoAneg9  = mean( Vel20_AoAneg9{2:end, 2:end-1} , 1);
avg_Vel20_AoAneg6  = mean( Vel20_AoAneg6{2:end, 2:end-1} , 1);
avg_Vel20_AoAneg3  = mean( Vel20_AoAneg3{2:end, 2:end-1} , 1);
avg_Vel20_AoA0     = mean( Vel20_AoA0{2:end, 2:end-1} , 1);
avg_Vel20_AoApos3  = mean( Vel20_AoApos3{2:end, 2:end-1} , 1);
avg_Vel20_AoApos6  = mean( Vel20_AoApos6{2:end, 2:end-1} , 1);
avg_Vel20_AoApos9  = mean( Vel20_AoApos9{2:end, 2:end-1} , 1);
avg_Vel20_AoApos12 = mean( Vel20_AoApos12{2:end, 2:end-1} , 1);
avg_Vel20_AoApos15 = mean( Vel20_AoApos15{2:end, 2:end-1} , 1);
avg_Vel20_AoApos18 = mean( Vel20_AoApos18{2:end, 2:end-1} , 1);
avg_Vel20_AoApos21 = mean( Vel20_AoApos21{2:end, 2:end-1} , 1);
avg_Vel20_AoApos24 = mean( Vel20_AoApos24{2:end, 2:end-1} , 1);


avg_Vel30_AoApos0  = mean(Vel30_AoApos0{2:end, 2:end-1}, 1);
avg_Vel30_AoApos3  = mean(Vel30_AoApos3{2:end, 2:end-1}, 1);
avg_Vel30_AoApos6  = mean(Vel30_AoApos6{2:end, 2:end-1}, 1);
avg_Vel30_AoApos9  = mean(Vel30_AoApos9{2:end, 2:end-1}, 1);
avg_Vel30_AoApos12 = mean(Vel30_AoApos12{2:end, 2:end-1}, 1);
avg_Vel30_AoApos15 = mean(Vel30_AoApos15{2:end, 2:end-1}, 1);
avg_Vel30_AoApos18 = mean(Vel30_AoApos18{2:end, 2:end-1}, 1);
avg_Vel30_AoApos21 = mean(Vel30_AoApos21{2:end, 2:end-1}, 1);
avg_Vel30_AoApos24 = mean(Vel30_AoApos24{2:end, 2:end-1}, 1);


% Coalesce the data into one table (matrix form)
Headers = Vel20_AoAneg12.Properties.VariableNames;
Vel20_data = [
    avg_Vel20_AoAneg12;
    avg_Vel20_AoAneg9;
    avg_Vel20_AoAneg6;
    avg_Vel20_AoAneg3;
    avg_Vel20_AoA0;
    avg_Vel20_AoApos3;
    avg_Vel20_AoApos6;
    avg_Vel20_AoApos9;
    avg_Vel20_AoApos12;
    avg_Vel20_AoApos15;
    avg_Vel20_AoApos18;
    avg_Vel20_AoApos21;
    avg_Vel20_AoApos24;
    ]

Headers = Vel30_AoApos12.Properties.VariableNames;
Vel30_data = [
    avg_Vel30_AoApos0;
    avg_Vel30_AoApos3;
    avg_Vel30_AoApos6;
    avg_Vel30_AoApos9;
    avg_Vel30_AoApos12;
    avg_Vel30_AoApos15;
    avg_Vel30_AoApos18;
    avg_Vel30_AoApos21;
    avg_Vel30_AoApos24;
    ]

%Solve for lift and drag force to plug into CL and CD equations
Lift_20 = cosd(Vel20_data(:,9)).*Vel20_data(:,1) + sind(Vel20_data(:,9)).*Vel20_data(:,2);
Drag_20 = sind(Vel20_data(:,9)).*Vel20_data(:,1) + cosd(Vel20_data(:,9)).*Vel20_data(:,2);

Lift_30 = cosd(Vel30_data(:,9)).*Vel30_data(:,1) + sind(Vel30_data(:,9)).*Vel30_data(:,2);
Drag_30 = sind(Vel30_data(:,9)).*Vel30_data(:,1) + cosd(Vel30_data(:,9)).*Vel30_data(:,2);


%Solve for CL and CD
CL_20 = Lift_20./(0.5.*Vel20_data(:,7).*Vel20_data(:,8).^2*A);
CD_20 = Drag_20./(0.5.*Vel20_data(:,7).*Vel20_data(:,8).^2*A);

CL_30 = Lift_30./(0.5.*Vel30_data(:,7).*Vel30_data(:,8).^2*A);
CD_30 = Drag_30./(0.5.*Vel30_data(:,7).*Vel30_data(:,8).^2*A);

[CL_20Max, CLmax20idx] = max(CL_20);
[CL_30Max, CLmax30idx]  = max(CL_30);
[CD_20Max, CDmax20idx] = max(CL_20);
[CD_30Max, CDmax30idx] = max(CL_30);

AoA_CL_20Max = 21.118;
AoA_CD_20Max = 21.118;
AoA_CL_30Max = 24.0547;
AoA_CD_30Max = 24.0547;

%Error Propagation of CL, CD, and RE. All values plus and minus
Force = 0.1; %N
TunnelSpeed = 0.4; %m/sec
AirDensity = 0.02;
AirViscosity = 0.01;
AirfoilDimensions = 0.001; %m
AoA = 0.2; %Degrees




%Plot CL, CD and CL/CD

ErrorMatrix = [];


figure(1)
plot(Vel20_data(:,9), CL_20, Vel30_data(:,9), CL_30)
xlabel('Angle of Attack (degrees)')
ylabel('C_L')
legend('20 m/s', '30 m/s')
grid on

figure(2)
plot(Vel20_data(:,9), CD_20, Vel30_data(:,9), CD_30)
xlabel('Angle of Attack (degrees)')
ylabel('C_D')
% errorbar(Vel20_data(:,9), CD_20)
legend('20 m/s', '30 m/s')
grid on

figure(3)
plot(Vel20_data(:,9), CL_20./CD_20, Vel30_data(:,9), CL_30./CD_30)
xlabel('Angle of Attack (degrees)')
ylabel('CL/CD')
legend('20 m/s', '30 m/s')
grid on