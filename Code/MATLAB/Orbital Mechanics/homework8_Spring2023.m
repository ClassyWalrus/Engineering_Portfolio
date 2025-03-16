close all;


% find places you need to replace
% by looking for "replace" variable
replace = 0;

% struct to save information
pathToFolder = 'C:\Users\davee\Desktop\Classes 2022\Orbital Mech\';
firstName = 'Colton'; % replace with your first and last name!
lastName = 'Davies'; % replace
HomeworkTitle = 'Homework8';
saveFilename = [HomeworkTitle '_' firstName '_' lastName '.mat'];
savePathAndFilename = [pathToFolder saveFilename];


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



%% part 1a 
aHohmann = (1/2)*(rHohmannPerigee+rHohmannApogee)
eccHohmann = (aHohmann-rHohmannPerigee)/aHohmann
%Depart from parking orbit about earth



% delta v is at periapsis where the velocity and delta v will  be
% tangential


% find the vel at apoapsis and periapsis of hohmann

% velocity magnitude at periapsis (departing from Earth parking orbit)
vHohmannPerigee = sqrt(2*((muEarth/rHohmannPerigee)-(muEarth/(aHohmann^2))))

% at apoapsis (arriving at lunar orbit)
vHohmannApoapsis = sqrt(2*((muEarth/rHohmannApogee)-(muEarth/(2*aHohmann))))

h8_p1a_struct = {aHohmann,eccHohmann,vHohmannPerigee,vHohmannApoapsis};
%% part 1b

% find the circular velcoity magnitude along this orbit (earth parking orbit)
v_park_Earth = sqrt(muEarth/rEarthParkingOrbit)

% delta v1 (departure from Earth parking orbit
delta_v_departEarthParkOrbit = abs(v_park_Earth-vHohmannPerigee)

h8_p1b_struct = {v_park_Earth,delta_v_departEarthParkOrbit};

%% part 1c
% find the circular velcoity along this orbit (lunar parking orbit)
v_park_Moon = sqrt(muMoon/rLunarParkingOrbit)

% what is the velocity of the moon with respect to earth
vMoon = sqrt(muEarth/aMoon)

% The spacecraft is slower than the moon

% Speed in moon-centered view at incoming asymptote (all velocities are
% tangential)

v_inf_m_Moon = abs(vHohmannApoapsis - vMoon)

% energy equation
v_hyperbola_atMoonParkOrbit = sqrt(v_inf_m_Moon^2+(2*muMoon/rLunarParkingOrbit))

% TO find delta v2, the deltav to transfer from the hyperbola at periapsis
% to the circular parking orbit about the moon

delta_v_hyperbolaToMoonParkOrbit =  v_hyperbola_atMoonParkOrbit - v_park_Moon

h8_p1c_struct = {v_park_Moon,vMoon,v_inf_m_Moon,v_hyperbola_atMoonParkOrbit,delta_v_hyperbolaToMoonParkOrbit};
%% part 1d
% find the total delta v
delta_v_total = delta_v_departEarthParkOrbit + delta_v_hyperbolaToMoonParkOrbit


% find the TOF
Period  = 2*pi*sqrt(aHohmann^3/muEarth)
TOF = Period/2
TOFdays = TOF/(86400)

% find period of moon
PeriodMoon = 2*pi*sqrt(aMoon^3/muEarth)
moonDegreesMoved = (vMoon/aMoon)*TOF


h8_p1d_struct = {delta_v_total,Period,TOF,TOFdays,PeriodMoon,moonDegreesMoved};

%% part 1e
% 
% Plotting
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
ypos_moon = rs_moon.*sin(thetas)
plot(xpos_moon,ypos_moon,'x','Color','y'); % 
legend('Earth','Earth parking orbit','Moon orbit');       
axis equal;


eccHohmann = 1 - rHohmannPerigee/aHohmann

rs = thetas*0;
for ri = 1:length(thetas)
    rs(ri) = aHohmann*(1-eccHohmann^2)/(1+eccHohmann*cos(thetas(ri)));
end

plot(rs(1:end*.5).*cos(thetas(1:end*.5)), rs(1:end*.5).*sin(thetas(1:end*.5)),'k-'); % plot ellipse
plot(rs(end*.5+1:end).*cos(thetas(end*.5+1:end)), rs(end*.5+1:end).*sin(thetas(end*.5+1:end)),'k-.'); % plot ellipse
axis equal;

legend('Earth','Earth parking orbit','Moon orbit', 'Hohmann Transfer','unused traj.'); 

%% Problem 2

% part a is vector diagram (can draw)
% part b is also a vector diagram

%% part c
%%

% use energy equation to find semi - major axis of new hyperbola

abs_a_new_Moon = muMoon/v_inf_m_Moon^2


% find new eccentricity 
ecc_new_Moon = rLunarParkingOrbit/abs_a_new_Moon+1

% The triangle formed by v_inf_m_Moon v_inf_p_Moon and delta_v_eq is an
% isosceles triangle (legs of hyperbola the same)

flyByAngle = 2*asin(1/ecc_new_Moon);


% now use the fact that we have an isosceles triangle to find delta_v_eq
delta_v_eq = 2*abs(v_inf_m_Moon)*sin(flyByAngle/2)


v_plus = sqrt(vMoon^2+v_inf_m_Moon^2-(2*vMoon*v_inf_m_Moon*cos(flyByAngle)))
r_plus = rHohmannApogee;

gammaNew = asin(v_inf_m_Moon*sin(flyByAngle)/v_plus)

energyNew = v_plus^2/2 - muEarth/aMoon

aNew = muEarth/(2*energyNew)

% h = rv*cos(gamma) = sqrt(mu p) -> (rv*cos(gamma))^2 /mu = p
pNew = (aMoon*v_plus*cos(gammaNew))^2/muEarth

eNew = sqrt(1+pNew/aNew)

r_periapsis = aNew*(energyNew-1)

thetaNew = acosd((pNew/aMoon-1)/eNew)

h8_p2c_struct = {abs_a_new_Moon,ecc_new_Moon,flyByAngle,delta_v_eq,v_plus,r_plus,gammaNew,energyNew,aNew...
    pNew,eNew,r_periapsis,thetaNew};

%% part 2d
% plot ellipse 2
thetas = ((thetaNew):1:50+thetaNew)*pi/180;
rs_new = thetas*0;
for ri = 1:length(thetas)
    rs_new(ri) =pNew/(1+eNew*cos(thetas(ri)));
end

% change in periapsis
delta_omega = -(thetaNew-180);
plot(rs_new.*cos(thetas + delta_omega*pi/180), rs_new.*sin(thetas+ delta_omega*pi/180),'-','Color','r'); % plot ellipse


legend('Earth','Earth parking orbit','Moon orbit', 'Hohmann Transfer','unused traj.', 'traj. after lunar flyby');



save(savePathAndFilename,"h8_p1a_struct","h8_p1b_struct","h8_p1c_struct","h8_p1d_struct","h8_p2c_struct");