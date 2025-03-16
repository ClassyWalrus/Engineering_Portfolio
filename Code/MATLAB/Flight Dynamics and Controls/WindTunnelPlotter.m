load("CRT_data_2024_3_28_16_50.mat")

run      =   data_matrix(:,1);
aoa      =   data_matrix(:,2);
beta     =   data_matrix(:,3);
elevator =   data_matrix(:,4);
rudder   =   data_matrix(:,5);
aileron  =   data_matrix(:,6);
density  =   data_matrix(:,7);
speed    =   data_matrix(:,8);
normal_F = - data_matrix(:,9);
axial_F  =   data_matrix(:,13);
normal_M =   data_matrix(:,17);

LE = .9;
TE = .9850;
c  = .2129;
S  = .5*c*(LE+TE);
xoverc = 0.25;

Q  = .5*density.*speed.^2;
L  = axial_F.*sind(aoa)-normal_F.*cosd(aoa);
D  = axial_F.*cosd(aoa)-normal_F.*sind(aoa);

CL = L./(Q.*S);
CD = D./(Q.*S);
CM = normal_M./(Q*S*c);

figure
scatter(aoa,CL,'filled')
grid on
grid minor
title('Coefficient of Lift for Varying Angle of Attack','FontSize',18)
ylabel('C_L','FontSize',16)
xlabel('\alpha (\circ)','FontSize',16)

% figure
% scatter(CL,CM,"filled")
% grid on
% grid minor
% title('Coefficient of Moment VS Coefficient of Lift','FontSize',18)
% xlabel('C_L','FontSize',16)
% ylabel('C_M ','FontSize',16)
% 
% figure
% scatter(CL,CM,"filled")
% grid on
% grid minor
% title('Coefficient of Moment VS Coefficient of Lift','FontSize',18)
% xlabel('C_L','FontSize',16)
% ylabel('C_M ','FontSize',16)

figure
scatter(aoa,CM,"filled")
grid on
grid minor
title('Coefficient of Moment for Varying Angle of Attack','FontSize',18)
xlabel('AoA','FontSize',16)
ylabel('C_M ','FontSize',16)

figure
scatter(CL,CD,"filled")
grid on
grid minor
title('Coefficient of Moment VS Coefficient of Lift','FontSize',18)
xlabel('C_L','FontSize',16)
ylabel('C_M ','FontSize',16)

