close all; clear all;
% find places you need to replace
% by looking for "replace" variable
replace = 0;

% struct to save information
pathToFolder = 'C:\Users\davee\Desktop\Classes 2022\Orbital Mech\';
firstName = 'Colton'; % replace with your first and last name!
lastName = 'Davies'; % replace
Title = 'midterm1';
saveFilename = [Title '_' firstName '_' lastName '.mat'];
savePathAndFilename = [pathToFolder saveFilename];

% constant values
R_Earth = 6378.137; % km
mu_Earth = 398600.4418; % km^3/s % note, mass of space vehicle is negligible

%% Problem 3
% values related to the question
altitude = 220; % km
r_earth = 6378.137; % km 
r_circular = r_earth + altitude; % 6578 km

% 3ai
% equation found in Textbook Eq 2.90
v_perigee = sqrt(2*mu_Earth/(r_circular)) % km/s

% 3aii
% equation found in Textbook Eq 2.63
v_circular = sqrt(mu_Earth/(r_circular)) % km/s

% 3aiii
% equation found in Textbook Eq 2.92 
% zero
zero = v_perigee/v_circular - sqrt(2)

% 3aiv
% equation found in Textbook Eq 2.90 - 2.93
dV_required = v_perigee-v_circular

h_m1_struct3a = {v_perigee,v_circular,zero,dV_required};

% 3b i.-v.
% equation found in Textbook Eq 2.90
v_i =  sqrt(2*mu_Earth/(2*r_earth)) % km/s
v_ii = sqrt(2*mu_Earth/(20*r_earth)) % km/s
v_iii = sqrt(2*mu_Earth/(30*r_earth)) % km/s
v_iv = sqrt(2*mu_Earth/(60*r_earth)) % km/s
v_v = 0 % km/s

h_m1_struct3b = {v_i,v_ii,v_iii,v_iv,v_v};

% setup for 3c (calculate r at perigee)
rperi_220 = (r_earth+220); % km
rperi_3000 = (r_earth+3000); % km
rperi_5000 = (r_earth+5000); % km

% setup for 3c (calculate semi latus rectum)
p220 = 2*rperi_220;
p3000 = 2*rperi_3000;
p5000 = 2*rperi_5000;

% setup for 3c (calculate r at 90 degree true anomaly)
r220_90d = p220;
r3000_90d = p3000;
r5000_90d = p5000;

% setup for 3c (calculate h)
h220 = sqrt(r220_90d*mu_Earth);
h3000 = sqrt(r3000_90d*mu_Earth);
h5000 = sqrt(r5000_90d*mu_Earth);

% setup for 3c (calculate v at 90 degree true anomaly)
v220 = sqrt(2*mu_Earth/(r220_90d)) %  km/s
v3000 = sqrt(2*mu_Earth/(r3000_90d)) % 
v5000 = sqrt(2*mu_Earth/(r5000_90d)) % 

% 3ci 
% equation found in Textbook Eq 2.93
gamma220 = 90/2
% 3cii
% equation found in Textbook Eq 2.93
gamma3000 = 90/2
% 3ciii
% equation found in Textbook Eq 2.93
gamma5000 = 90/2

h_m1_struct3c = {gamma220,gamma3000,gamma5000};

% PLOTTING %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
thetas = -125:1:125; % theta values

% r values at the thetas above
rs220 = p220./(1+cosd(thetas));
rs3000 = p3000./(1+cosd(thetas));
rs5000 = p5000./(1+cosd(thetas));

color220 = 'k'; color3000 = 'b'; color5000 = 'r';
figure;
% Parabola 220
plot(rs220.*cosd(thetas), rs220.*sind(thetas), 'Color',color220); hold on;
plot(rs3000.*cosd(thetas), rs3000.*sind(thetas), 'Color',color3000); hold on;
plot(rs5000.*cosd(thetas), rs5000.*sind(thetas), 'Color',color5000); hold on;
legend('parabola perigee alt.= 220 km','parabola perigee alt.= 3000 km','parabola perigee alt.= 5000 km','AutoUpdate','off')
% Position Vector at theta = 90
theta = 90;
plot(0, r220_90d,'x','Color',color220);
text(0, r220_90d, '\leftarrow satellite @ \theta =90\circ' );
plot([0 -.65*r220_90d*cosd(theta)],[0 .5*r220_90d*sind(theta)], '-.', 'Color',color220); % r
text(-.65*r220_90d*cosd(theta), .5*r220_90d*sind(theta), '\leftarrow r' );

% Velocity Vector
ecc = 1;
v_vec220 = mu_Earth/h220*[(-sind(theta))   (ecc+cosd(theta))]; 
quiver(r220_90d*cosd(theta),r220_90d*sind(theta),v_vec220(1)*rperi_220/norm(v_vec220),v_vec220(2)*rperi_220/norm(v_vec220) ,'Color', color220); % plot velocity
text(r220_90d*cosd(theta)+v_vec220(1)*.5*rperi_220/norm(v_vec220),r220_90d*sind(theta)+v_vec220(2)*.5*rperi_220/norm(v_vec220), '  \leftarrow vdir');

% plot the local horizon 
q1a = quiver(r220_90d*cosd(theta),r220_90d*sind(theta),rperi_220*sind(theta),-rperi_220*cosd(theta) ,'Color',color220); % plot local horizon
q1a.ShowArrowHead = 'off';
q1a.LineStyle = '-.';


% Directrix
xline(2*rperi_220,'-','Color',color220);
text(2*rperi_220, 0, '\leftarrow directrix');

xlabel('$$\hat{p} (km)$$','Interpreter','Latex');
ylabel('$$\hat{q} (km)$$','Interpreter','Latex');


%%%  Plotting for 3000 km altitude
% NEED TO DO ALL THE PLOTTING FOR 3000 km
% Position Vector at theta = 90
theta = 90;
plot(0, r3000_90d,'x','Color',color3000);
text(0, r3000_90d, '\leftarrow satellite @ \theta =90\circ' );
plot([0 -.65*r3000_90d*cosd(theta)],[0 .5*r3000_90d*sind(theta)], '-.', 'Color',color3000); % r
text(-.65*r3000_90d*cosd(theta), .5*r3000_90d*sind(theta), '\leftarrow r' );

% Velocity Vector
ecc = 1;
v_vec3000 = mu_Earth/h3000*[(-sind(theta))   (ecc+cosd(theta))]; 
quiver(r3000_90d*cosd(theta),r3000_90d*sind(theta),v_vec3000(1)*rperi_3000/norm(v_vec3000),v_vec3000(2)*rperi_3000/norm(v_vec3000) ,'Color', color3000); % plot velocity
text(r3000_90d*cosd(theta)+v_vec3000(1)*.5*rperi_3000/norm(v_vec3000),r3000_90d*sind(theta)+v_vec3000(2)*.5*rperi_3000/norm(v_vec3000), '  \leftarrow vdir');

% plot the local horizon 
q1a = quiver(r3000_90d*cosd(theta),r3000_90d*sind(theta),rperi_3000*sind(theta),-rperi_3000*cosd(theta) ,'Color',color3000); % plot local horizon
q1a.ShowArrowHead = 'off';
q1a.LineStyle = '-.';


% Directrix
xline(2*rperi_3000,'-','Color',color3000);
text(2*rperi_3000, 0, '\leftarrow directrix');


%%%  Plotting for 5000 km altitude
% NEED TO DO ALL THE PLOTTING FOR 5000 km
% Position Vector at theta = 90
theta = 90;
plot(0, r5000_90d,'x','Color',color5000);
text(0, r5000_90d, '\leftarrow satellite @ \theta =90\circ' );
plot([0 -.65*r5000_90d*cosd(theta)],[0 .5*r5000_90d*sind(theta)], '-.', 'Color',color5000); % r
text(-.65*r5000_90d*cosd(theta), .5*r5000_90d*sind(theta), '\leftarrow r' );

% Velocity Vector
ecc = 1;
v_vec5000 = mu_Earth/h5000*[(-sind(theta))   (ecc+cosd(theta))]; 
quiver(r5000_90d*cosd(theta),r5000_90d*sind(theta),v_vec5000(1)*rperi_5000/norm(v_vec5000),v_vec5000(2)*rperi_5000/norm(v_vec5000) ,'Color', color5000); % plot velocity
text(r5000_90d*cosd(theta)+v_vec5000(1)*.5*rperi_5000/norm(v_vec5000),r5000_90d*sind(theta)+v_vec5000(2)*.5*rperi_5000/norm(v_vec5000), '  \leftarrow vdir');

% plot the local horizon 
q1a = quiver(r5000_90d*cosd(theta),r5000_90d*sind(theta),rperi_5000*sind(theta),-rperi_5000*cosd(theta) ,'Color',color5000); % plot local horizon
q1a.ShowArrowHead = 'off';
q1a.LineStyle = '-.';


% Directrix
xline(2*rperi_5000,'-','Color',color5000);
text(2*rperi_5000, 0, '\leftarrow directrix');

%% Problem 4
eccentricities = [.05 .1 .15 .2 .25 .3]*pi;
thetas_E90 = [acos((eccentricities(1,1)-cosd(90))/(eccentricities(1,1)*cosd(90)-1)) acos((eccentricities(1,2)-cosd(90))/(eccentricities(1,2)*cosd(90)-1)) acos((eccentricities(1,3)-cosd(90))/(eccentricities(1,3)*cosd(90)-1)) acos((eccentricities(1,4)-cosd(90))/(eccentricities(1,4)*cosd(90)-1)) acos((eccentricities(1,5)-cosd(90))/(eccentricities(1,5)*cosd(90)-1)) acos((eccentricities(1,6)-cosd(90))/(eccentricities(1,6)*cosd(90)-1))];
figure; plot(-cos(thetas_E90), eccentricities);
xlabel('$$ -Cos(theta) (radians)$$', 'Interpreter', 'Latex');
ylabel('$$ Eccentricities (unitless)$$', 'Interpreter', 'Latex');
title('Eccentricities vs. -cos(theta)');

h_m1_struct4 = {thetas_E90};
save(savePathAndFilename,"h_m1_struct3a","h_m1_struct3b","h_m1_struct3c","h_m1_struct4");

           