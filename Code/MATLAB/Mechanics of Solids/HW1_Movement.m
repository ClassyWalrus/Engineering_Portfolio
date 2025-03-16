%%Colton Davies
%%Mechanics of Solids
%%HW1
clear all;
close all;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%
%%1. Initialize spatial and temporal grid
%%Setup spatial geometry
w = 2; % [m] sample width
L = 1; % [m] sample length
c = -1;
a = 2;
b = 1;
d = 1;
v = 2;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%
%%2. Initialize test parameters
simTime = 2; % [s] total simulation time
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%
%%3. Setup spatial grid
nrZ = 20; % number of gridpoints along z axis
dZ = L/(nrZ-1); % spacing between gridpoints
z = 0:dZ:L; % z position of each gridpoint
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%
%%4. Set up temporal grid
dT = 0.5; % [s] time steps for algorithm
nrT = round(simTime/dT); % number of time steps in simulation
Time = 0:dT:simTime;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%
%%Apply displacement equation at each time instance nrt
% for i = 1:length(Time)
% for j = 1:nrZ
% u(i,j) = (z(1,j)-c.*Time(1,i)).^3; %input equation for 2A

% end
% end
% for i = 1:length(Time)
% for j = 1:nrZ
% u(i,j) = Time(1,i)/2*a*(z(1,j)/Time(1,i))^2+b*Time(1,i)+c;
% % input equation for 2B
% end
% end
%
% u = u(sum(isnan(u),2)==0,:); %check u for NaN row
for i = 1:length(Time)
 for j = 1:nrZ
 u(i,j) = z(1,j)*(1+v*Time(1,i)/L); %input equation for 2A
 end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%
%%Create plot visual
% width of sample
y = w.* ones(1, length(z));
% plot size of specimen
%first plot location without strain
for s = 1:height(u)
 subplot(height(u),1,s)
 plot(u(s, :), y(1, :));
 sgtitle('y(x, t)')
% plot vertical lines
line_coordinates = [u(s,:)' u(s,:)' zeros(size(z')) y'];
line_coordinates = mat2cell(line_coordinates, ones(numel(z),1), [2 2])';
line_coordinates(3,:) = {'k'};
plot(line_coordinates{:})
line([u(s,1),u(s,nrZ)],[w,w])
ylim([0 1.5*w])
xlim([0 max(u(:,nrZ))])
end