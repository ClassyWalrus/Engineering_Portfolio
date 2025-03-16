Zeta = (sqrt(2)*0.035)/(2*0.75);
wn = sqrt(2)*(32.2/745);
wd = wn*sqrt(1-zeta^2);


A1 = -44.01;
A2 = -0.083057;
A2_Gamma = -0.0022751;


t = 0:0.1:300;
U = (exp(-wn.*Zeta.*t)).*(A1.*cos(wd.*t)+A_2.*sin(wd.*t));
Gamma =(exp(-wn.*Zeta.*t)).*(A2_Gamma.*sin(wd.*t));


figure
plot(t,U)
ylabel('Velocity (ft/s)')
xlabel('Time (seconds)')
title('Velocity Perturbation Response U(t)')


figure
plot(t,Gamma)
ylabel('Flight Path Angle (radians)')
xlabel('Time (seconds)')
title('Flight Path Angle Perturbation Response ')