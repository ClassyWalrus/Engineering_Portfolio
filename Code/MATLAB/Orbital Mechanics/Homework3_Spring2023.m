% Homework 3 Problem 2

% find places you need to replace
% by looking for "replace" variable
replace = 0;

% struct to save information
pathToFolder = 'C:\Users\davie378\Documents\MATLAB\';
firstName = 'Colton'; % replace with your first and last name!
lastName = 'Davies'; % replace
HomeworkTitle = 'Homework3';
saveFilename = [HomeworkTitle '_' firstName '_' lastName '.mat'];
savePathAndFilename = [pathToFolder saveFilename];


% Problem Values
v_vec =  [.5 3 0]; % km/s
mass_SC = 420; % kg
alt_SC = 8530; % km

% constants from Curtis
G = (6.6742 * 10^-11)*10^-9;% Newtons kg-2 m2
mass_Earth = 5.974*10^24; %kg
radius_Earth = 6378; % km

% assign variables related to two-body problem
m1 = mass_Earth; %kg
m2 = mass_SC; % kg
mu = G*(m1+m2);
r = 8550 + radius_Earth; % km

% figure out time derivative variables from velocity vector
r_dot = v_vec(1); 
theta_dot = v_vec(2)/r;

r_vec = [r, 0, 0];

% part a
% Equation from Lecture 5 Slide 4 ((m1*m2)/(m1+m2))*(r x r_dot);
C3 = ((m1*m2)/(m1+m2))*cross([r,0,0],v_vec);

disp(['The answer to Part A is [' num2str(C3) '] km^2*kg/s'])

% part b
% Equation from Lecture 6 Slide 3 (r x r_dot)

h = cross(r_vec,v_vec);

disp(['The answer to Part B is ' num2str(h(3)) ' km^2/s'])

% part c
% Equations from Lecture 6 Slide 4
T = 0.5*(m1*m2)/(m1+m2)*dot(v_vec,v_vec); % Kinetic Energy
U = G*m1*m2/r; % Potential Energy
C4 = T-U; % System Energy

disp(['The answer to Part C is ' num2str(C4) ' kg*km^2/s^2'])

% part d
% Equation from Lecture 6 Slide 4
specificEnergy = 0.5*dot(v_vec,v_vec)-mu/r;

disp(['The answer to Part D is ' num2str(specificEnergy) ' km^2/s^2'])

% part e
% Equation from Lecture 6 Slide 4
dAdt = h/2;

disp(['The answer to Part E is ' num2str(dAdt(3)) ' km^2/s'])

% part f
% Equation Lecture 5 Slide 11
p = h(3)^2/mu;

disp(['The answer to Part F is ' num2str(p) ' km'])

% part g
% Equation from Lecture 5 Slide 11
e = sqrt(specificEnergy*2*h(3)^2/mu^2+1);

disp(['The answer to Part G is ' num2str(e) ' '])

% part h
% Equation from Lecture 5 Slide 11
a = p/(1-e^2);

disp(['The answer to Part F is ' num2str(a) ' km'])

% part i
% Equation from Lecture 5 Slide 4
Period = 2*pi*sqrt(a^3/mu);

disp(['The answer to Part I is ' num2str(Period) ' s'])

% part j
% Equation from 5 Slide 11
theta = acos((p/r-1)/e);

disp(['The answer to Part J is ' num2str(theta) ' radians'])

% part k
% Equation from 
rp = r*cos(theta);
rq = r*sin(theta);

disp(['The answer to Part K is ' num2str(rp) ' km and ' num2str(rq) ' km'])

% Plotting the orbit
thetas = (1:1:360)*pi/180; % radians
rs = zeros(length(thetas),1); % initialize
for ri = 1:length(thetas)
    rs(ri) = p/(1+e*cos(thetas(ri)));
end


% %%%%%%%%%%%%%%%%%%%%%%% PLOTTING %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure(1); plot(rs.*cos(thetas)', rs.*sin(thetas)','k-'); % plot ellipse
hold on;
plot([0 0],[0 p], 'm-.'); % plot semi latus rectum
text(0,p/2, '\leftarrow p' );
plot(-a*e, 0, 'kx'); % plot the center
text(-a*e, 0, '\leftarrow C' );

plot(0,0,'go'); % attracting body
text(0,0, '\leftarrow attracting body' );
text(0,p/2, '\leftarrow p' );
rth = p/(1+e*cos(theta));
plot(rth*cos(theta), rth*sin(theta),'mx');
text(rth*cos(theta), rth*sin(theta), '\leftarrow satellite' );

title('Orbit Plane Plot');
xlabel('x (km)'); ylabel('y (km)');
axis equal;
% %%%%%%%%%%%%%%%%%%%%%%% END PLOTTING %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


h3struct = {C3, h, T, U, C4, specificEnergy,...
            dAdt, p, e, a, Period,theta,rp, rq};
save(savePathAndFilename,"h3struct");

