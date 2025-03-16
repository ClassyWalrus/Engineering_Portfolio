% Aerodynamic numbers
Cd = 0.04;
s = 0.017;
m = 0.003;
g = 9.81;
Cl = 0.2205;
p = 1.225;
tspan = [0, 6];

% Initial Conditions
x_a = [4; 0; 2; 0];
x_b = [11; 0; 2; 0];

% Solve ODEs using anonymous function to pass extra parameters
[a_t, a_x] = ode45(@(t,x) airplaneEqs(t, x, Cd, s, m, g, Cl, p), tspan, x_a);
[b_t, b_x] = ode45(@(t,x) airplaneEqs(t, x, Cd, s, m, g, Cl, p), tspan, x_b);

% Plot results
figure;
plot(a_t, a_x(:,1), b_t, b_x(:,1))
title('Airspeed (U) vs Time');
legend('Airspeed:4m/s','Airspeed:11m/s');

figure;
plot(a_t, a_x(:,2), b_t, b_x(:,2));
title('Flight Path Angle (Γ) vs Time');
legend('Airspeed:4m/s','Airspeed:11m/s');

figure;
plot(a_t, a_x(:,3), b_t, b_x(:,3));
title('Altitude (h) vs Time');
legend('Airspeed:4m/s','Airspeed:11m/s');

figure;
plot(a_t, a_x(:,4), b_t, b_x(:,4));
title('Distance (d) vs Time');
legend('Airspeed:4m/s','Airspeed:11m/s');

function dxdt = airplaneEqs(t, x, Cd, s, m, g, Cl, p)
    U = x(1);
    G = x(2);
    U_dot = -(Cd*p*s*U^2)/(2*m) - g*sin(G);
    G_dot = (Cl*p*s*U)/(2*m) - (g/U)*cos(G);
    h_dot = U*sin(G);
    d_dot = U*cos(G);
    dxdt = [U_dot; G_dot; h_dot; d_dot];
end
