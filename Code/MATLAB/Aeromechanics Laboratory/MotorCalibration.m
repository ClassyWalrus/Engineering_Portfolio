% This was a record of some of the live calculations and checks done on the
% SYSTEMS 1 Lab (Thrust Coefficient Section)

% Simeon Shaffar
% Connor Morley
% Colton Davies

clc; close all; clear;

%% Setting up the data
drone_mass = 1.175; % kg

data_table = readtable("Log_2024-10-24_092042.csv");

input_signal_table = data_table{:, "ESCSignal__s_"}; % Digital signal

thrust_kgf = data_table{:, "Thrust_kgf_"}; % kgf
thrust_N = thrust_kgf * 9.81; % N

rotor_rpm = data_table{:, "MotorOpticalSpeed_RPM_"}; % RPM
rotor_rad_per_s = rotor_rpm * 2*pi/60; % rad/s


%% Interpreting the data
kt = 5E-6; % N/(rad/s)^2
curve_fit = @(w) kt*w^2;
[P,S] = polyfit(rotor_rad_per_s.^2,thrust_N,1);
% R_sq = 1 - (S.normr/norm(E_Vals - mean(E_Vals)))^2
f1 = polyval(P,rotor_rad_per_s.^2);
disp(S)

figure
hold on 
grid on
plot(rotor_rad_per_s, thrust_N, '.b', MarkerSize=10)
fplot(curve_fit, [0, 1000], '--r', LineWidth=1.5)

title("Thrust vs rotor speed of drone propeller")
xlabel("\omega (rad/s)")
ylabel("Thrust (N)")
legend("Measured Data", "Curve fit (T = 5E-6*\omega^2)")

figure(2)
hold on 
grid on
plot(rotor_rad_per_s.^2, thrust_N, '.b', MarkerSize=10)
plot(rotor_rad_per_s.^2, f1, '--r', LineWidth=1.5)

title("Thrust vs rotor speed of drone propeller")
xlabel("\omega^2 (rad/s)")
ylabel("Thrust (N)")
legend("Measured Data", "Curve fit (T = 5.14E-6*\omega^2)")