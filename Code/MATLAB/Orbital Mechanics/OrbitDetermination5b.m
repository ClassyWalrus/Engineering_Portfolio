load("OrbitDetermineFinalProject.mat");

%plot3(Xtot(:,1),Xtot(:,2),Xtot(:,3),'r.');
CoeFromSV = zeros(345:7);
for i = 1:1:345
    R_Vec = [Xtot(i,1) Xtot(i,2) Xtot(i,3)];
    V_Vec = [Xtot(i,4) Xtot(i,5) Xtot(i,6)];
    CoeFromSV(i,:) = coe_from_sv(R_Vec, V_Vec, 398600);
end

CoeSemi = CoeFromSV(:,7).*(1-CoeFromSV(:,2).^2);
Mu_Earth = 398600;
J2 = 1.08263e-3;

GPS_TIME_Sized = GPS_TIME(1:344,2);
timeInSeconds = GPS_TIME(:,2) - GPS_TIME(1,2);
timeInSecondsDiff = diff(GPS_TIME);
R = CoeFromSV(:,7).*(1-CoeFromSV(:,2).^2)./(1+CoeFromSV(:,2).*cosd(CoeFromSV(:,6)));

%Eccentricites
figure
subplot(3,2,1);
plot(timeInSeconds,CoeFromSV(:,2))
xlabel('GPS Time (s)'); ylabel('Eccentricity');


%RAAN
%figure(2);
subplot(3,2,2);
plot(timeInSeconds,CoeFromSV(:,3))
xlabel('GPS Time (s)'); ylabel('RAAN (degrees)');


%Inclination
%figure(3);
subplot(3,2,3);
plot(timeInSeconds,CoeFromSV(:,4))
xlabel('GPS Time (s)'); ylabel('Inclination (degrees)');


%Argument of Perigee
%figure(4);
subplot(3,2,4);
plot(timeInSeconds,CoeFromSV(:,5))
xlabel('GPS Time (s)'); ylabel('Argument of Perigee (degrees)');


%Semimajor axis
%figure(5);
subplot(3,2,5);
plot(timeInSeconds,CoeFromSV(:,7))
xlabel('GPS Time (s)'); ylabel('Semimajor axis (km)');


%Semilatus rectum
%figure(6);
subplot(3,2,6);
plot(timeInSeconds, CoeSemi)
xlabel('GPS Time (s)'); ylabel('Semilatus rectum (km)');

sgtitle('Orbital Parameters v.s. Time(s)')


%Time rate of change for RAAN using J2
figure(2);
dRAAN_J2 = -((3/2)*(sqrt(Mu_Earth)*J2*R)/(((1-CoeFromSV(:,2).^2).^2).*CoeFromSV(:,2).^(7/2)))*(cosd(CoeFromSV(:,4)));

plot(timeInSeconds,dRAAN_J2, 'r')
xlabel('GPS Time (s)'); ylabel('Derivative of RAAN');
title('Time rate of change for RAAN using J2')

 
% % Time rate of change for Argument of Perigee using J2
figure(3);
dPerigee_J2 = -((3/2)*(sqrt(Mu_Earth)*J2*R)/(((1-CoeFromSV(:,2).^2).^2).*CoeFromSV(:,2).^(7/2)))*((5/2)*(sind(CoeFromSV(:,4)).^2)-2);

plot(timeInSeconds,dPerigee_J2, 'k')
xlabel('GPS Time (s)'); ylabel('Argument of Perigee Derivative');
title('Time rate of change for Argument of Perigee using J2')

%Time rate of change for RAAN using diff
figure(4);
dRAAN=diff(CoeFromSV(:,3))./diff(GPS_TIME(:,2));

plot(timeInSeconds(2:1:345),dRAAN, 'r')
xlabel('GPS Time (s)'); ylabel('Derivative of RAAN');
title('Time rate of change for RAAN using diff')

%Time rate of change for Argument of Perigee using diff
figure(5);
dPerigee = diff(CoeFromSV(:,5))./diff(GPS_TIME(:,2));

plot(timeInSeconds(2:1:345),dPerigee,'k')
xlabel('GPS Time (s)'); ylabel('Argument of Perigee Derivative');
title('Time rate of change for Argument of Perigee using diff')


