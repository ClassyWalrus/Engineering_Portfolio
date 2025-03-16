close all; clear all;

% find places you need to replace
% by looking for "replace" variable
replace = 0;

% struct to save information
pathToFolder = 'C:\Users\davee\Desktop\Classes 2022\Orbital Mech';
firstName = 'Colton'; % replace with your first and last name!
lastName = 'Davies'; % replace
HomeworkTitle = 'Homework7';
saveFilename = [HomeworkTitle '_' firstName '_' lastName '.mat'];
savePathAndFilename = [pathToFolder saveFilename];


%% Problem 1 - Option 1
Mu_Earth = 398600.4418; % km^2/s^2
R_Earth =  6378.137; % km
ecc = 0.65;
a =  5.8*R_Earth;
% Lecture 8 Slide 8
p = a*(1-ecc^2); % semi latus rectum

% using the relationship between semi-major axis and period
wait_time = (2*pi)*sqrt(a^3/(Mu_Earth))/60/60; % half a period (in hours)

% using geometry of an ellipse
%Lecture 19 Slide 8
r1 = a*(1+ecc); %km

% using the specific energy to find velocity before maneuver
v1_minus = sqrt(2*(((-Mu_Earth)/(2*a))+(Mu_Earth/r1))); %km/s

% using the specific energy to find velocity after maneuver
v1_plus = sqrt(2*((Mu_Earth/a)-(Mu_Earth/(2*a))));
% Take the difference to find delta V
deltaV1 = v1_plus - v1_minus


% plot ellipse
thetas = (-180:1:180)*pi/180; % thetas in radians!
rs = thetas*0;
for ri = 1:length(thetas)
    
    theta_ri = thetas(ri);
    % Use conic equation to calculate r values
    rs(ri) = p/(1+ecc*cos(theta_ri));
end

figure; plot(rs.*cos(thetas), rs.*sin(thetas),'k-'); % plot ellipse
hold on;

% Plot r at specific true anomaly (theta) of maneuver for option 1 (apoapsis) in
% radians
theta = pi;

% use conic equation to calculate r at  theta
rth = p/(1+ecc*cos(theta));

plot(rth*cos(theta), rth*sin(theta),'mx');
text(rth*cos(theta), rth*sin(theta), '\leftarrow satellite @ t_a' );
plot([0 rth*cos(theta)],[0 rth*sin(theta)], 'b-.'); % r
text(.5*rth*cos(theta), .5*rth*sin(theta), '\leftarrow r_a' );


% plot the velocity direction using quiver at r_a
v_veca = [0 -v1_minus];
quiver(rth*cos(theta),rth*sin(theta),v_veca(1)*.5*a/norm(v_veca),v_veca(2)*.5*a/norm(v_veca) ,'k'); % plot velocity
text(rth*cos(theta)+v_veca(1)*.25*a/norm(v_veca),rth*sin(theta)+v_veca(2)*.25*a/norm(v_veca), '  \leftarrow vdir');


% Plot the circle
thetas = (-180:1:180)*pi/180;
rs = thetas*0;
for ri = 1:length(thetas)
    rs(ri) = r1;
end

plot(rs.*cos(thetas), rs.*sin(thetas),'k-'); % plot ellipse
hold on;

% Plot r at specific theta (theta a)
rth = r1;
plot(rth*cos(theta), rth*sin(theta),'mx');
text(rth*cos(theta), rth*sin(theta), '\leftarrow satellite @ t_a' );
plot([0 rth*cos(theta)],[0 rth*sin(theta)], 'b-.'); % r
text(.5*rth*cos(theta), .5*rth*sin(theta), '\leftarrow r_a' );

% plot the velocity direction using quiver at r_a
v_veca = [0 -v1_plus]; % in perifocal frame
quiver(rth*cos(theta),rth*sin(theta),v_veca(1)*.5*a/norm(v_veca),v_veca(2)*.5*a/norm(v_veca) ,'k'); % plot velocity
text(rth*cos(theta)+v_veca(1)*.25*a/norm(v_veca),rth*sin(theta)+v_veca(2)*.25*a/norm(v_veca), '  \leftarrow vdir');



title('Orbits in perifocal frame');
xlabel('p (km)'); ylabel('q (km)');
axis equal

h7_p1_opt1_struct = {p, wait_time, r1, v1_minus, v1_plus, deltaV1, theta, rth};

%% Problem 1  - Option 2
thetad = 125;
r1 = p/(1+ecc*cosd(thetad)); % conic equation

E = acosd((a-r1)/(a*ecc));
% E is between 0 and 180 because theta is between 0 and 180

% time since periapsis
t1_mtp = ((E-sind(E))/(sqrt(Mu_Earth/a^3)))/60/60;

% we can use v_1m = rdot *rhat + r*thetadot*thetahat
v1m_mag = sqrt(2*Mu_Earth/r1 - Mu_Earth/a)

h = sqrt(Mu_Earth*p);

% components of velocity in perifocal lecture 19 slide 6
v_r = ((Mu_Earth*ecc)/h)*sind(thetad);
v_theta = h/r1;

% flight path angle lecture 19 slide 9
gamma = acosd((sqrt(Mu_Earth*p))/(r1*v1m_mag));


% To circularize, we know that the new velocity must be in the 
% theta hat direction (gamma = 0)
v1p_mag = sqrt(Mu_Earth/r1)

gammaNew = acosd((sqrt(Mu_Earth*p))/(r1*v1p_mag));

% Lecture 19 slide 10
vN = v1p_mag; va = v1m_mag;
deltaV = sqrt(vN^2+va^2-2*vN*va*cosd(gammaNew));

% turning angle
alpha = -(180- acosd( (vN^2 - deltaV^2 - va^2) / (-2*deltaV*va)))



% plot ellipse
thetas = (-180:1:180)*pi/180;
rs = thetas*0;
for ri = 1:length(thetas)
    rs(ri) = p/(1+ecc*cos(thetas(ri)));
end

figure; plot(rs.*cos(thetas), rs.*sin(thetas),'k-'); % plot ellipse
hold on;
% Plot r at specific theta (theta a)
theta = pi/180*thetad;
rth = p/(1+ecc*cos(theta));
plot(rth*cos(theta), rth*sin(theta),'mx');
text(rth*cos(theta), rth*sin(theta), '\leftarrow satellite @ t_a' );
plot([0 rth*cos(theta)],[0 rth*sin(theta)], 'b-.'); % r
text(.5*rth*cos(theta), .5*rth*sin(theta), '\leftarrow r_a' );

% plot the velocity direction using quiver at r_a
v_veca = [-sind(thetad) (ecc+cosd(thetad))]*Mu_Earth/h;
quiver(rth*cos(theta),rth*sin(theta),v_veca(1)*.5*a/norm(v_veca),v_veca(2)*.5*a/norm(v_veca) ,'k'); % plot velocity
text(rth*cos(theta)+v_veca(1)*.25*a/norm(v_veca),rth*sin(theta)+v_veca(2)*.25*a/norm(v_veca), '  \leftarrow vdir');


% Plot the circle
thetas = (-180:1:180)*pi/180;
rs = thetas*0;
for ri = 1:length(thetas)
    rs(ri) = r1;
end

plot(rs.*cos(thetas), rs.*sin(thetas),'k-'); % plot ellipse
hold on;
% Plot r at specific theta (theta a)
theta = pi/180*thetad;
rth = r1;
plot(rth*cos(theta), rth*sin(theta),'mx');
text(rth*cos(theta), rth*sin(theta), '\leftarrow satellite @ t_a' );
plot([0 rth*cos(theta)],[0 rth*sin(theta)], 'b-.'); % r
text(.5*rth*cos(theta), .5*rth*sin(theta), '\leftarrow r_a' );

% plot the velocity direction using quiver at r_a
v_veca = [-sind(thetad) cosd(thetad)]*Mu_Earth/h;
quiver(rth*cos(theta),rth*sin(theta),v_veca(1)*.5*a/norm(v_veca),v_veca(2)*.5*a/norm(v_veca) ,'k'); % plot velocity
text(rth*cos(theta)+v_veca(1)*.25*a/norm(v_veca),rth*sin(theta)+v_veca(2)*.25*a/norm(v_veca), '  \leftarrow vdir');



title('Orbits in perifocal frame');
xlabel('p (km)'); ylabel('q (km)');
axis equal



h7_p1_opt2_struct = { r1, E, t1_mtp, h, v_r, v_theta, gamma,v1p_mag, gammaNew,deltaV}

%% Problem 2
rc = 5.7*R_Earth;
vc = 3.4; % km/s 
gammac = -10;

% Lecture 19 , Slide 10
ecc = sqrt( (rc*vc^2/Mu_Earth -1)^2*(cosd(gammac))^2 + (sind(gammac))^2)
thetac = atand(((rc*vc^2)/(Mu_Earth))*cosd(gammac)*sind(gammac))/(((rc*vc^2)/(Mu_Earth))*cos(gammac)^2-1);

hc = rc*vc*cosd(gammac);
p = hc^2/Mu_Earth;
a = p/(1-ecc^2)

Period = 2*pi*sqrt(a^3/Mu_Earth)/3600  

specificEnergy = -Mu_Earth/(2*a)

% periapsis distance
rp = a*(1-ecc)

% apoapsis distance
ra = a*(1+ecc)

% Note that since the flight path angle is negative, the E must be in
% negative portion of orbit
Ec = acos((a-rc)/(a*ecc));
Ec = 2*pi-Ec
Ec_deg = Ec*180/pi

tc_m_tp = sqrt(a^3/Mu_Earth)*(Ec - ecc*sin(Ec))/3600

h7_p2_parta_struct = { thetac, hc, p, a, rp, ra};
%% 2b

rm = p
Em = acos((a-rm)/a*ecc); % eccentric anomaly in radians
Em_deg = Em*180/pi % eccentric anomaly in degrees

tm_m_tp = sqrt(a^3/Mu_Earth)*(Em - ecc*sin(Em))/3600
% Find time until maneuver
tm_m_tc = Period - tc_m_tp +tm_m_tp 

% figure out velocity at maneuver
vm = sqrt(2*Mu_Earth/rm - Mu_Earth/a)

% flight path angle
gamma = acosd(sqrt(Mu_Earth*p)/(rm*vm));

deltaV = 1; % kms
alpha  = 30;

vN = sqrt(deltaV^2 + vm^2 -2*deltaV*vm*cosd(180-alpha))



deltaGamma = acosd( (deltaV^2 - vm^2 - vN^2)/(2*vm*vN))
% since delta gamma is turning away from Earth
deltaGamma = 180- deltaGamma

% new flight path angle
gammaN = gamma + deltaGamma

% mow calculate new values after maneuver
eccN = sqrt( (rm*vN^2/Mu_Earth -1)^2*(cosd(gammaN))^2 + (sind(gammaN))^2);
thetaN = atand(((rm*vN^2)/(Mu_Earth))*cosd(gammaN)*sind(gammaN))/(((rm*vN^2)/(Mu_Earth))*cosd(gammaN)^2-1);

% new magnitude of specific angular momentum 
% new semi-latus rectum, new semi-major axis
hN = rm*vN*cosd(gammaN);
pN = hN^2/Mu_Earth;
aN = pN/(1-eccN^2);

PeriodN = 2*pi*sqrt(aN^3/Mu_Earth)/3600

specificEnergyN = -Mu_Earth/(2*aN)

% new periapsis distance
rpN = aN*(1-eccN)

% new apoapsis distance
raN = aN*(1+eccN);

% new eccentric anomaly
EN = acos((aN-rm)/(aN*eccN))
tn_m_tp = sqrt(aN^3/Mu_Earth)*(EN-eccN*sin(EN))/3600

% figure out what happened to our periapsis angle (how much did it change)
delta_omega = 90 - thetaN


% plot ellipse 1
thetas = (-180:1:180)*pi/180;
rs = thetas*0;
for ri = 1:length(thetas)
    rs(ri) = p/(1+ecc*cos(thetas(ri)));
end

figure; plot(rs.*cos(thetas), rs.*sin(thetas),'k-'); % plot ellipse
hold on;


% plot ellipse 2
thetas = (-180:1:180)*pi/180;
rsN = thetas*0;
for ri = 1:length(thetas)
    rsN(ri) = pN/(1+eccN*cos(thetas(ri)));
end

plot(rsN.*cos(thetas + delta_omega*pi/180), rsN.*sin(thetas+ delta_omega*pi/180),'r-'); % plot ellipse
legend('Old','New');
xlabel('Old p hat (km)');
ylabel('Old q hat (km)');
axis equal

h7_p2_partb_struct = { Em, gamma, eccN, thetaN, hN, pN, aN,raN};

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

save(savePathAndFilename,"h7_p1_opt1_struct","h7_p1_opt2_struct","h7_p2_parta_struct","h7_p2_partb_struct");