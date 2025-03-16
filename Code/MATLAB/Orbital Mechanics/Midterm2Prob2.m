replace = 0;
aMoon = 384400; % km
r_earth = 6378; % km
r_moon = 1737; %km
muEarth = 398600; % km^3/s^2
muMoon = 4902.5; % km^3/s^2

parkOrbitEarthAlt = 210; % km
parkOrbitMoonAlt = 160; % km
rLunarParkingOrbit = r_moon + parkOrbitMoonAlt; % km
rEarthParkingOrbit = r_earth + parkOrbitEarthAlt; %km
rHohmannApogee = aMoon - rLunarParkingOrbit %km
rHohmannPerigee = rEarthParkingOrbit % km
aHohmann = (1/2)*(rHohmannPerigee+rHohmannApogee)



%knowns from HW8
rOld = 382503; %km
P = 6588; %km

%Part B) i-iv
ecc = ((P/rOld)-1)/cosd(175);
a = P/(1-ecc^2);
rApogee = P/(1+ecc*cosd(180)) % km
v = sqrt(2*((muEarth/rApogee)-(muEarth/(2*a))))
Gamma = acosd(sqrt(muEarth*P)/(rApogee*v))

%Part D)
vEarth = sqrt(muEarth/rApogee)
vInf = sqrt(vEarth^2+v^2-2*vEarth*v*cosd(Gamma));

%Part E)
%Using geometry
FlyByAngle = 2*asin(1/ecc);

%Part F)
rp = rLunarParkingOrbit;

%Part G) i-vi
SpecEnergy = (vInf^2)/2;
aNew = replace;
PNew = replace;
OrbitalRadius = replace;
NewTrueAnomaly = replace;
DeltaOmega = replace;



%PLOTTING
figure; 
plot(0,0,'bx'); legend('Earth');
hold on;
xlabel('p direction (km)');
ylabel('q direction(km)');

thetas = (1:1:360)*pi/180;

rs_earthParking = rEarthParkingOrbit*ones(size(thetas));
xpos_sat = rs_earthParking.*cos(thetas );
ypos_sat = rs_earthParking.*sin(thetas);
plot(xpos_sat,ypos_sat,'.','Color','b'); % plot parking orbit
legend('Earth','Earth parking orbit');       
axis equal;

rs_moon = aMoon*ones(size(thetas));
xpos_moon = rs_moon.*cos(thetas);
ypos_moon = rs_moon.*sin(thetas);
plot(xpos_moon,ypos_moon,'x','Color','y'); % 
legend('Earth','Earth parking orbit','Moon orbit');       
axis equal;

eccHohmann = 1 - rHohmannPerigee/aHohmann

rs = thetas*0;
for ri = 1:length(thetas)
    rs(ri) = aHohmann*(1-eccHohmann^2)/(1+eccHohmann*cos(thetas(ri)));
end

plot(rs(1:end*.5).*cos(thetas(1:end*.5)), rs(1:end*.5).*sin(thetas(1:end*.5)),'k-'); % plot ellipse
axis equal;

rs = thetas*0;
for ri = 1:length(thetas)
    rs(ri) = a*(1-ecc^2)/(1+ecc*cos(thetas(ri)));
end

plot(rs(1:end*.5).*cos(thetas(1:end*.5)), rs(1:end*.5).*sin(thetas(1:end*.5)),'g'); % plot ellipse
plot(rs(end*.5+1:end).*cos(thetas(end*.5+1:end)), rs(end*.5+1:end).*sin(thetas(end*.5+1:end)), 'Color','#EDB120'); % plot ellipse
axis equal;

legend('Earth','Earth parking orbit','Moon orbit', 'Hohmann Transfer','unused traj.');
