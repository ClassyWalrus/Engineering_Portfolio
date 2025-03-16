clc; clear; close all;

%% A Construction
u_o = 824.2;
h = 33000;
rho = 0.000795;
W = 230000;
g = -32.2;
theta_o = 0;
b = 142.3;
c = 23;
m = 7143;
S = 2600;
V_To = 824.2;

Ixz = 45000;
Ix = 3770000;
Iz = 7130000;
star_coef = 1/(1-(Ixz^2)/(Ix*Iz));

C_y_beta = -0.7277;
C_l_beta = -0.16732;
C_l_p = -0.516;
C_l_r = 0.147;
C_n_beta = 0.15471;
C_n_p = -0.0107;
C_n_r = -0.190;
C_n_d_a = -0.003701;
C_n_d_r = -0.08337;
C_l_d_a = -0.07965;
C_l_d_r = 0.021086;
C_y_d_a = 0.18651;

Y_v = (rho*S*V_To/(2*m))*C_y_beta;
Y_p = 0;
Y_r = 0;

L_beta = (rho*S*V_To^2*b/(2*Ix))*C_l_beta;
L_p = (rho*S*V_To*b^2/(4*Ix))*C_l_p;
L_r = (rho*S*V_To*b^2/(4*Ix))*C_l_r;

N_beta = (rho*S*V_To^2*b/(2*Iz))*C_n_beta;
N_p = (rho*S*V_To*b^2/(4*Iz))*C_n_p;
N_r = (rho*S*V_To*b^2/(4*Iz))*C_n_r;

L_beta_star = L_beta*star_coef;
L_p_star = L_p*star_coef;
L_r_star = L_r*star_coef;

N_beta_star = N_beta*star_coef;
N_p_star = N_p*star_coef;
N_r_star = N_r*star_coef;

%% B construction

Y_d_r = ((rho*S*V_To^2)/(2*m))*C_y_d_a;

L_d_a = ((rho*S*V_To^2*b)/(2*Ix))*C_l_d_a;
L_d_r = ((rho*S*V_To^2*b)/(2*Ix))*C_l_d_r;

N_d_a = ((rho*S*V_To^2*b)/(2*Iz))*C_n_d_a;
N_d_r = ((rho*S*V_To^2*b)/(2*Iz))*C_n_d_r;

L_d_a_star = L_d_a*star_coef;
L_d_r_star = L_d_r*star_coef;

N_d_a_star = N_d_a*star_coef;
N_d_r_star = N_d_r*star_coef;

%% A + B Complete

A = [Y_v    Y_p    Y_r-u_o    g*cos(theta_o);
     (L_beta_star + Ixz/Ix*N_beta_star)/u_o    L_p_star + Ixz/Ix*N_p_star    L_r_star + Ixz/Ix*N_r_star    0;
     (N_beta_star + Ixz/Iz*L_beta_star)/u_o    N_p_star + Ixz/Iz*L_p_star    N_r_star + Ixz/Iz*L_r_star    0;
     0 1 0 0];

B = [0    Y_d_r;
     L_d_a_star + Ixz/Ix*N_d_a_star    L_d_r_star + Ixz/Ix*N_d_r_star;
     N_d_a_star + Ixz/Iz*L_d_a_star    N_d_r_star + Ixz/Iz*L_d_r_star;
     0    0];

%% EqMotion Initial Conditions
v_o = 16;    %Side Slip Velocity
p_o = 0.3;   %Roll Rate
r_o = 0;   %Yaw Rate
phi_o = 0.1; %Roll Angle

x0 = [v_o ; p_o ; r_o ; phi_o];
%% Open Poles

OpenPoles = eig(A);
% Open_Combined = A-B*OpenPoles;
%% Closed Poles
zeta = 0.99;
omega_n = 0.499;
omega_3 = 0.499;
omega_4 = 0.45;

ClosedPoles = ones(1,4);
ClosedPoles(1) = -zeta*omega_n+i*omega_n*sqrt(1-zeta^2);
ClosedPoles(2) = -zeta*omega_n-i*omega_n*sqrt(1-zeta^2);
ClosedPoles(3) = -omega_3;
ClosedPoles(4) = -omega_4;

K = place(A,B,ClosedPoles);

Closed_Combined = A-B*K;


%% Graphing

%Open Params
figure
grid on
[t,x] = ode45(@(t,x) odefunc2(t,x, A), [0:0.001:30],x0);
subplot(4,1,1)
plot(t',x(:,1)');
grid on
title('Side Slip Velocity')
xlabel('Time(s)')
ylabel('v(t) (ft/s)')

subplot(4,1,2)
plot(t',x(:,2)')
grid on
title('Roll Rate')
xlabel('Time(s)')
ylabel('p(t) (rad/s)')

subplot(4,1,3)
plot(t',x(:,3)');
grid on
title('Yaw Rate')
xlabel('Time(s)')
ylabel('r(t) (rad/s)')

subplot(4,1,4)
plot(t',x(:,4)');
grid on
title('Roll Angle')
xlabel('Time(s)')
ylabel('\phi(t) (rad)')

sgtitle('Open Loop')

%Closed Params
figure
[t,x] = ode45(@(t,x) odefunc(t,x, Closed_Combined), [0:0.001:30],x0);
subplot(4,1,1)
plot(t,x(:,1));
grid on
title('Side Slip Velocity')
xlabel('Time(s)')
ylabel('v(t) (ft/s)')

subplot(4,1,2)
plot(t,x(:,2));
grid on
title('Roll Rate')
xlabel('Time(s)')
ylabel('p(t) (rad/s)')

subplot(4,1,3)
plot(t,x(:,3));
grid on
title('Yaw Rate')
xlabel('Time(s)')
ylabel('r(t) (rad/s)')

subplot(4,1,4)
plot(t,x(:,4));
grid on
title('Roll Angle')
xlabel('Time(s)')
ylabel('\phi(t) (rad)')

sgtitle('Closed Loop')

%Open and Closed Poles
figure
hold on
grid on
plot(real(OpenPoles),imag(OpenPoles), 'ko')
plot(real(ClosedPoles),imag(ClosedPoles), 'ro')
legend('Open','Closed')
title('Complex Conjugate Poles')
xlabel('Real Values')
ylabel('Imaginary Values')
text(real(ClosedPoles(1))+0.05, imag(ClosedPoles(1)), 'Dutch')
text(real(ClosedPoles(2))+0.05, imag(ClosedPoles(2)), 'Dutch')
text(real(ClosedPoles(3))+0.05, imag(ClosedPoles(3)), 'Roll')
text(real(ClosedPoles(4))+0.05, imag(ClosedPoles(4)), 'Spiral')

text(real(OpenPoles(1))+0.05, imag(OpenPoles(1)), 'Dutch')
text(real(OpenPoles(2))+0.05, imag(OpenPoles(2)), 'Dutch')
text(real(OpenPoles(3))+0.05, imag(OpenPoles(3)), 'Roll')
text(real(OpenPoles(4))+0.05, imag(OpenPoles(4)), 'Spiral')


%Control Surface Deflection
Deflection = [];
N = size(t);
for i = 1:N
    Deflection(i,:) = -K*x(i,:)'*(180/pi);
end

figure
hold on

subplot(2,1,1)
plot(t,Deflection(:,1));
grid on
title('Aileron Deflection')
xlabel('Time (s)')
ylabel('\delta_\alpha (deg)')

subplot(2,1,2)
plot(t,Deflection(:,2));
grid on
title('Rudder Deflection')
xlabel('Time (s)')
ylabel('\delta_r (deg)')

sgtitle('Time History of Control Surface Deflection')

function y = odefunc(t,x,Closed_Combined)
    y = Closed_Combined*x;
end

function y = odefunc2(t,x,A)
    y = A*x;
end
