close all; clear all;

% find places you need to replace
% by looking for "replace" variable
replace = 0;

% struct to save information
pathToFolder = 'C:\Users\davee\Desktop\Classes 2022\Orbital Mech\Homework4_Spring2023.m';
firstName = 'Colton'; % replace with your first and last name!
lastName = 'Davies'; % replace
HomeworkTitle = 'Homework4';
saveFilename = [HomeworkTitle '_' firstName '_' lastName '.mat'];
savePathAndFilename = [pathToFolder saveFilename];


R_Earth = 6378.137; % km
mu_Earth = 398600.4418; % km^3/s % note, mass of space vehicle is negligible


% values related to the question
r_perigee = 4.5*R_Earth; % km
r_apogee = 22*R_Earth; % km
theta_deg = 115; % degrees
theta_rad = theta_deg*pi/180; % radians



% part a1
% equation from Lecture 6 Slide 13
a = (r_perigee+r_apogee)/2;


% part a2
% equation from Lecture 6 Slide 13
ecc = 1 - (r_perigee/a);

% part a3
% equation from Lecture 6 Slide 13?
p = r_perigee*(1+ecc);

% part a4
% equation from?
b = sqrt(a^2-(a*ecc)^2);

% part a5
% equation from Textbook, Eq 2.70?
specificEnergy = -mu_Earth/(2*a);

% part a6
% equation from Textbook, Eq 2.61?
h = sqrt(mu_Earth*a*(1-ecc^2));

% part a7
% equation from Textbook, Eq 2.73?
P_second = (2*pi)/(sqrt(mu_Earth))*a^(3/2);
P_hour = P_second/(60*60);

% part a8
% equation from Textbook, Eq 2.78?
r_a = sqrt(r_apogee*r_perigee);

% part a9
% equation from?
v_a = replace;

% part a10
% equation from Textbook, Eq 2.42?
gamma_deg_a = tan((ecc*sin(theta_rad))/(1+ecc*cos(theta_rad)));
gamma_rad_a = gamma_deg_a*pi/180;

% part a11
% equation from Textbook, Eq 2.53?
v_circ = sqrt(mu_Earth/((r_apogee+r_perigee)/2));



%% b

% time that has passed since periapsis
deltaT_b = 21*60*60; % time pass in seconds

% part b1
% equation from?
E_rad_a = replace;
E_deg_a = 180/pi*E_rad_a;

% part b2 (make sure you are using radians!!)
% equation from?
M_a = replace;

% part b3
% equation from?
n = replace;

% part b4
% equation from?
tsp = replace;


% part b5
% equation from?
v_vec_a = [replace  replace];

% part b6
% equation ?
M_b = replace;
M_deg_b = 180*M_b/pi

% part b7
% equation from?
E_b =kepler_E(replace, replace);
E_deg_b = E_b*180/pi

% part b8
% equation from?
r_b = replace;

% part b9
% equation from?
v_b = replace;

% part b10
% equation from?
gamma_b = replace;
if (E_b > pi)
    gamma_b = -1*gamma_b;
end
gamma_deg_b = gamma_b*180/pi

% part b11
% equation from?
theta_b = replace;
theta_deg_b = theta_b*180/pi;

% part b12
% equation from?
v_vec_b = [replace  replace];

% part b13
% equation from?
r_p = r_b*cos(theta_b)
r_q = replace;



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% PLOTTING %%%%%%%%%%%%%%%%%%%%%%
thetas = (1:1:360)*pi/180;
rs = thetas*0;
for ri = 1:length(thetas)
    rs(ri) = p/(1+ecc*cos(thetas(ri)));
end

% plot the ellipse, semi latus rectum and ellipse center, and attracting
figure; plot(rs.*cos(thetas), rs.*sin(thetas),'k-'); % plot ellipse
hold on;
plot([0 0],[0 p], 'm-.'); % plot semi latus rectum
text(0,p/2, '\leftarrow p' );
plot(-a*ecc, 0, 'kx'); % plot the center
text(-a*ecc, 0, '\leftarrow C' );
plot(0,0,'go'); % attracting body
text(0,0, '\leftarrow attracting body' );

% Plot r at specific theta (theta a)
rth = p/(1+ecc*cos(theta_rad));
plot(rth*cos(theta_rad), rth*sin(theta_rad),'mx');
text(rth*cos(theta_rad), rth*sin(theta_rad), '\leftarrow satellite @ t_a' );
plot([0 rth*cos(theta_rad)],[0 rth*sin(theta_rad)], 'b-.'); % r
text(.5*rth*cos(theta_rad), .5*rth*sin(theta_rad), '\leftarrow r_a' );

% plot the velocity direction using quiver at r_a
quiver(rth*cos(theta_rad),rth*sin(theta_rad),v_vec_a(1)*.5*a/norm(v_vec_a),v_vec_a(2)*.5*a/norm(v_vec_a) ,'k'); % plot velocity
text(rth*cos(theta_rad)+v_vec_a(1)*.25*a/norm(v_vec_a),rth*sin(theta_rad)+v_vec_a(2)*.25*a/norm(v_vec_a), '  \leftarrow vdir');

% plot the local horizon 
q1a = quiver(rth*cos(theta_rad),rth*sin(theta_rad),-a*sin(theta_rad),a*cos(theta_rad) ,'k'); % plot local horizon
q1a.ShowArrowHead = 'off';
q1a.LineStyle = '-.';

q2a = quiver(rth*cos(theta_rad),rth*sin(theta_rad),a*sin(theta_rad),-a*cos(theta_rad) ,'k'); %  plot local horizon
q2a.ShowArrowHead = 'off';
q2a.LineStyle = '-.';
text(rth*cos(theta_rad)+.5*a*sin(theta_rad),rth*sin(theta_rad)-.5*a*cos(theta_rad), '\leftarrow local horizon');

% Plot r at specific theta (theta b)
% NEED TO DO PLOTTING AT theta b
replace;
replace;
replace;

title('Orbit Plane Plot');
xlabel('x (km)'); ylabel('y (km)');
axis equal;

%%%%%%%%%%%%%%%%%%%%%% END PLOTTING %%%%%%%%%%%%%%%%%%%%%%%%%%%%
h4struct = {a, ecc, p, b, specificEnergy, h, P_hour, r_a, v_a, gamma_deg_a, v_circ, ...
deltaT_b, E_deg_a, M_a, n, tsp, v_vec_a, M_deg_b, E_deg_b, r_b, v_b, gamma_deg_b,...
theta_deg_b, v_vec_b, r_p, r_q};
save(savePathAndFilename,"h4struct");