load("CRT_data_2024_3_27_14_39.mat")

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
W = 7.48;

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

figure
scatter(CL,CM,"filled")
grid on
grid minor
title('Coefficient of Moment VS Coefficient of Lift','FontSize',18)
xlabel('C_L','FontSize',16)
ylabel('C_M ','FontSize',16)

figure
scatter(aoa,CM,"filled")
grid on
grid minor
title('Coefficient of Moment for Varying Angle of Attack','FontSize',18)
xlabel('\alpha (\circ)','FontSize',16)
ylabel('C_M ','FontSize',16)

%%Problem 2 Solutions
ClMax = max(CL)
StaticMargin = -((0.075168--0.0751031)/(-0.0701657-0.819761))
CL_AlphaSlope = (0.787817--0.898177)/(10.05--10.047)
PitchStiffness = (0.430829-0.261908)/(-10.059-8.013)
X_NPmeters = (xoverc+StaticMargin)*c
X_NPfrac = X_NPmeters/c
X_NPinches = (xoverc+StaticMargin)*c*39.37
ElevatorPower =(0.387046--0.280546)/(-36.5-38)


CL_MinV = W/(0.5*mean(density)*(8^2)*S)
CL_MaxV = W/(0.5*mean(density)*(25^2)*S)
CM_Naut = 0.05


ForwardLimit = X_NPmeters + (c/CL_MinV)*(CM_Naut + ElevatorPower*-18)
AftLimit = X_NPmeters + (c/CL_MaxV)*(CM_Naut + ElevatorPower*18)

