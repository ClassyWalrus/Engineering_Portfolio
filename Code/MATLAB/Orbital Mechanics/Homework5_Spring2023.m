close all; clear all;
% find places you need to replace
% by looking for "replace" variable
replace = 0;

% struct to save information
pathToFolder = 'C:\Users\davee\Desktop\Classes 2022\Orbital Mech\';
firstName = 'Colton'; % replace with your first and last name!
lastName = 'Davies'; % replace
HomeworkTitle = 'Homework5';
saveFilename = [HomeworkTitle '_' firstName '_' lastName '.mat'];
savePathAndFilename = [pathToFolder saveFilename];

% constant values
R_Earth = 6378.137; % km
mu_Earth = 398600.4418; % km^3/s % note, mass of space vehicle is negligible
mu_Moon =  4903; % km^3/s^2 %note, mass of space vehicle is negligible
RMoon =  1737; % km

% values related to the question
theta = -65; % degree
absa = 7100; % km
a = -7100; % km
passageAlt = 800; % km



% part a1
% equation from geometry
r_perilune = RMoon+passageAlt;


% part a2
% equation Lecture 7 Slide 6
v_inf = sqrt(mu_Moon/abs(a));

% part a3
% equation Lecture 7 Slide 6
specificEnergy = v_inf^2/2;

% part a4
% equation Lecture 7 Slide 6
ecc =  1 - (r_perilune/a)

% part a5
% equation Lecture 7 Page 7
v_perilune =  sqrt(2*mu_Moon/RMoon);

% part a6
% equation Textbook Eq 2.96
h =  sqrt((r_perilune*mu_Moon)*(1 + ecc*cosd(theta)));

% part a7
% equation Lecture 7 Slide 7
delta =  2*asind(1/ecc);

% part a8
% equation Textbook Eq 2.106
b =  a*sqrt(ecc^2-1);

% part a9
% equation Lecture 7 Slide 7
p =  a*(1-ecc^2);

% part a10
% equation Lecture 10 Slide 9
r =  p/(1+ecc*cosd(theta));

% part a11
% equation Textbook Eq 2.49 and Eq 2.31
v = sqrt(((mu_Moon/h)*ecc*sind(theta))^2+(h/r)^2);

% part a12
% equation 
gamma =  atand(((mu_Moon/h)*ecc*sind(theta))/(h/r));

% part a13
% equation Lecture 10 Slide 10
F =  acosh(((r/absa+1)/ecc));

% part a14
% equation Lecture 7 Slide 7
thetaInf = 90+delta/2;

% part a15
% equation Lecture 7 Slide 6
Center_point_dist = absa+r_perilune;


%%
% part b
% equation 
thetas = (-120:1:120)*pi/180;
rs = thetas*0;
for ri = 1:length(thetas)
    rs(ri) =  p./(1+ecc*cos(thetas(ri)));
end


figure; plot(rs.*cos(thetas), rs.*sin(thetas),'k-'); % plot ellipse
hold on;

title('Orbit in perifocal frame');
xlabel('p (km)'); ylabel('q (km)');


theta = theta * pi/180;
% Plot r at specific theta (theta a)
r_ta =  p./(1+ecc*cos(theta));
plot(r_ta*cos(theta), r_ta*sin(theta),'mx');
text(r_ta*cos(theta), r_ta*sin(theta), '\leftarrow satellite @ t_a' );
plot([0 r_ta*cos(theta)],[0 r_ta*sin(theta)], 'b-.'); % r
text(-.65*r_ta*cos(theta), .5*r_ta*sin(theta), 'r_a\rightarrow ' );

% plot the velocity direction using quiver at r_a
v_veca = [(mu_Moon/h)*-sin(theta) (mu_Moon/h)*(ecc+cos(theta))];  
quiver(r_ta*cos(theta),r_ta*sin(theta),v_veca(1)*a/norm(v_veca),v_veca(2)*a/norm(v_veca) ,'k'); % plot velocity
text(r_ta*cos(theta)+v_veca(1)*.5*a/norm(v_veca),r_ta*sin(theta)+v_veca(2)*.5*a/norm(v_veca), '  \leftarrow vdir');

% plot the local horizon 
q1a = quiver(r_ta*cos(theta),r_ta*sin(theta),-a*sin(theta),a*cos(theta) ,'k'); % plot local horizon
q1a.ShowArrowHead = 'off';
q1a.LineStyle = '-.';

% plot the local horizon 
q1a = quiver(r_ta*cos(theta),r_ta*sin(theta),a*sin(theta),-a*cos(theta) ,'k'); % plot local horizon
q1a.ShowArrowHead = 'off';
q1a.LineStyle = '-.';
text(r_ta*cos(theta)+v_veca(1)*.5*a/norm(v_veca),r_ta*sin(theta)+v_veca(2)*.25*a/norm(v_veca), '  \leftarrow \gamma');
axis equal;


% plot the asymptotes
hold on; 
q1a = quiver(Center_point_dist,0,a*2*cosd(thetaInf), a*2*sind(thetaInf) ,'k'); 
q1a.ShowArrowHead = 'off';
q1a.LineStyle = '-.';

q1a = quiver(Center_point_dist,0,-a*2*cosd(thetaInf), -a*2*sind(thetaInf) ,'k'); % 
q1a.ShowArrowHead = 'off';
q1a.LineStyle = '-.';

q1a = quiver(Center_point_dist,0,a*2*cosd(thetaInf), -a*2*sind(thetaInf) ,'k'); 
q1a.ShowArrowHead = 'off';
q1a.LineStyle = '-.';

q1a = quiver(Center_point_dist,0,-a*2*cosd(thetaInf), a*2*sind(thetaInf) ,'k'); 
q1a.ShowArrowHead = 'off';
q1a.LineStyle = '-.';


q1a = quiver(Center_point_dist,0,0, -a*2 ,'k'); 
q1a.ShowArrowHead = 'off';
q1a.LineStyle = '-.';

q1a = quiver(Center_point_dist,0,0, a*2 ,'k'); 
q1a.ShowArrowHead = 'off';
q1a.LineStyle = '-.';

text(Center_point_dist-abs(a)*.25, abs(a)*.25, '\delta /2');

%% plotting the semi minor
halfDeltaFBA = (asin(1/ecc));
plot([0 b*cos(.5*pi-halfDeltaFBA)],[0 b*sin(.5*pi-halfDeltaFBA)],'k'); % plot local horizon
text( b*cos(.5*pi-halfDeltaFBA), b*sin(.5*pi-halfDeltaFBA),'b');


plot( b*cos(.5*pi-halfDeltaFBA),- b*sin(.5*pi-halfDeltaFBA),'rx'); % plot local horizon
text( b*cos(.5*pi-halfDeltaFBA), -b*sin(.5*pi-halfDeltaFBA),'aimpoint');
axis equal;
 

%% Problem 2
R_Earth = 6378.137; % km
a = 8*R_Earth; % km
theta_deg = 135; % deg
ecc = .45;
incl = 32; % deg
Omega = 40; % deg
omega = 65; % deg

r_ta = a*(1-ecc^2)/(1+ecc*cosd(theta_deg));
p = a*(1-ecc^2);


thetas = (1:1:360)*pi/180;
rs = thetas*0;
for ri = 1:length(thetas)
    rs(ri) = p/(1+ecc*cos(thetas(ri)));
end

%%%%%%%%%%%%%%%%%%%%%% PLOTTING %%%%%%%%%%%%%%%%%%%%%%%%%%%%
% plot the ellipse, semi latus rectum and ellipse center, and attracting
figure; plot(rs.*cos(thetas), rs.*sin(thetas),'k-'); % plot ellipse
hold on;
plot([0 0],[0 p], 'm-.'); % plot semi latus rectum
text(0,p/2, '\leftarrow p' );
plot(-a*ecc, 0, 'kx'); % plot the center
text(-a*ecc, 0, '\leftarrow C' );
plot(0,0,'go'); % attracting body
text(0,0, '\leftarrow attracting body' );
plot([0 r_ta*cosd(theta_deg)],[0 r_ta*sind(theta_deg)], 'b-.'); % r
text(.65*r_ta*cosd(theta_deg), .5*r_ta*sind(theta_deg), 'r_a' );
xlabel('p (km)'); ylabel('q (km)');
%%%%%%%%%%%%%%%%%%%%%% END PLOTTING %%%%%%%%%%%%%%%%%%%%%%%%%%%%

% rotation matrix
sO = sind(Omega); cO = cosd(Omega); si = sind(incl); ci = cosd(incl); so = sind(omega); co = cosd(omega);
rotmZYX = [-sO*ci*so+cO*co   -sO*ci*co-cO*so   sO*si;
           cO*ci*so+sO*co   cO*ci*co-sO*so   -cO*si;
           si*so   si*co   ci];

rs = zeros(length(thetas),3);
for ri = 1:length(thetas)
    r = p/(1+ecc*cos(thetas(ri)));
    rcos = r*cos(thetas(ri));
    rsin = r*sin(thetas(ri));
    rv =  rotmZYX*[rcos; rsin; 0];
    rs(ri,:) = rv;
end
rthi = rotmZYX*[r_ta*cosd(theta_deg); r_ta*sind(theta_deg); 0];
r_p = a*(1-ecc^2)/(1+ecc*cosd(0));
per_i = rotmZYX*[r_p*cosd(0); r_p*sind(0); 0];

%%%%%%%%%%%%%%%%%%%%%% PLOTTING %%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure;
plot3(rs(:,1), rs(:,2),rs(:,3),'k-'); % plot ellipse
hold on;

plot3([0 rthi(1)], [0 rthi(2)],[0 rthi(3)],'k-'); % plot r
plot3([0 per_i(1)], [0 per_i(2)],[0 per_i(3)],'b-'); % plot p
rectangle('Position',[-10*R_Earth -10*R_Earth 20*R_Earth 20*R_Earth], 'FaceColor',[.5 0 1 .5] );
view(110,25)

replace = '';
xlabel('Xhat (km)'); ylabel('Yhat (km)'); zlabel('Zhat (km)');

axis equal;
%%%%%%%%%%%%%%%%%%%%%% END PLOTTING %%%%%%%%%%%%%%%%%%%%%%%%%%%%

h5struct = {r_perilune,v_inf,specificEnergy,ecc,v_perilune,h,...
            delta, b,p,r,v, gamma, F};
save(savePathAndFilename,"h5struct");
           