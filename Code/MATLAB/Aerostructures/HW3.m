%%Aerospace Structures HW3
%%Colton Davies

%% Problem 2a

syms B;
L = 3;
E = 180e9;
rho = 3;
Ix = 5.375e-8;

A = [1 0 1 0;
     0 B 0 B;
     B^2*cosh(L*B) B^2*sinh(L*B) -B^2*cos(L*B) -B^2*sin(L*B);
     B^3*sinh(L*B) B^3*cosh(L*B) B^3*sin(L*B) -B^3*cos(L*B)]

detA = det(A)
f = @(B) B^6*(cos(3*B)^2 + cosh(3*B)^2 + 2*cos(3*B)*cosh(3*B) + sin(3*B)^2 - sinh(3*B)^2)

%found zeroes of the detA visually along graph
fplot(detA,[0,10])
xlabel('Beta')
ylabel('Determinant Values')
hold on
grid on

B1 = fzero(f,0.6);
B2 = fzero(f,1.5);
B3 = fzero(f,2.64);
B4 = fzero(f,3.65);
B5 = fzero(f,4.65);

B = [B1 B2 B3 B4 B5]

w = sqrt((E*Ix)/rho).*B.^2

%% Problem 2b
%For N=6
N6 = 6; N8 = 8; N10 = 10; N15 = 15; N30 = 30;
NodePos6 = (0:L/(N6-1):L)'; NodePos8 = (0:L/(N8-1):L)'; NodePos10 = (0:L/(N10-1):L)'; NodePos15 = (0:L/(N15-1):L)'; NodePos30 = (0:L/(N30-1):L)';
NE6 = N6-1; NE8 = N8-1; NE10 = N10-1; NE15 = N15-1; NE30 = N30-1;
ElmConnect6 = [1:1:N6-1; 2:1:N6]'; ElmConnect8 = [1:1:N8-1; 2:1:N8]'; ElmConnect10 = [1:1:N10-1; 2:1:N10]'; ElmConnect15 = [1:1:N15-1; 2:1:N15]'; ElmConnect30 = [1:1:N30-1; 2:1:N30]';
E = 180e9; E6 = ones(NE6,1)*E; E8 = ones(NE8,1)*E; E10 = ones(NE10,1)*E; E15 = ones(NE15,1)*E; E30 = ones(NE30,1)*E;
Ix = 5.375e-8; Ix6 = ones(NE6,1)*Ix; Ix8 = ones(NE8,1)*Ix; Ix10 = ones(NE10,1)*Ix; Ix15 = ones(NE15,1)*Ix; Ix30 = ones(NE30,1)*Ix;
rho = 3; rho6 = ones(NE6,1)*rho; rho8 = ones(NE8,1)*rho; rho10 = ones(NE10,1)*rho; rho15 = ones(NE15,1)*rho; rho30 = ones(NE30,1)*rho;
LEndNode = 1;
REndNode6 = 6; REndNode8 = 8; REndNode10 = 10; REndNode15 = 15; REndNode30 = 30;
Nmodes = 5;

%call gmeshbeammodes for all N values

[Freq6, RawModes6, M6, K6, InterpModes6, Zpos6] = gmesh_beam_modes(...
                                                      N6, NodePos6, ...
                                                      NE6, ElmConnect6, ...
                                                      E6, Ix6, rho6, ...
                                                      LEndNode, ...
                                                      REndNode6, ...
                                                      Nmodes);

[Freq8, RawModes8, M8, K8, InterpModes8, Zpos8] = gmesh_beam_modes(...
                                                      N8, NodePos8, ...
                                                      NE8, ElmConnect8, ...
                                                      E8, Ix8, rho8, ...
                                                      LEndNode, ...
                                                      REndNode8, ...
                                                      Nmodes);

[Freq10, RawModes10, M10, K10, InterpModes10, Zpos10] = gmesh_beam_modes(...
                                                      N10, NodePos10, ...
                                                      NE10, ElmConnect10, ...
                                                      E10, Ix10, rho10, ...
                                                      LEndNode, ...
                                                      REndNode10, ...
                                                      Nmodes);

[Freq15, RawModes15, M15, K15, InterpModes15, Zpos15] = gmesh_beam_modes(...
                                                      N15, NodePos15, ...
                                                      NE15, ElmConnect15, ...
                                                      E15, Ix15, rho15, ...
                                                      LEndNode, ...
                                                      REndNode15, ...
                                                      Nmodes);

[Freq30, RawModes30, M30, K30, InterpModes30, Zpos30] = gmesh_beam_modes(...
                                                      N30, NodePos30, ...
                                                      NE30, ElmConnect30, ...
                                                      E30, Ix30, rho30, ...
                                                      LEndNode, ...
                                                      REndNode30, ...
                                                      Nmodes);

%All frequencies for all Ns match hand calc.

%Plotting for first part of b
figure (2)
nodesnums = [6, 8, 10, 15, 30]
frequencymatrix = [Freq6, Freq8, Freq10, Freq15, Freq30]
hold on
plot(nodesnums, frequencymatrix(1,:),'k')
plot(nodesnums,frequencymatrix(2,:),'r')
plot(nodesnums,frequencymatrix(3,:),'g')
plot(nodesnums,frequencymatrix(4,:),'c')
plot(nodesnums,frequencymatrix(5,:),'m')
title('Frequencies vs Nodes')
legend('Freq1','Freq2','Freq3','Freq4','Freq5','location','northoutside')
xlabel('Number of Nodes')
ylabel('Vibration Frequency Value (rads/s)')
hold off

%Matlab has a weird bug with AMD graphics cards that doesn't allow a legend
%to color code.
%%Plotting for second part of b

figure (3)
hold on
plot(Zpos8,InterpModes8)
title('Normal Mode Shapes')
xlabel('Distance (m)')
ylabel('InterpModes8')
legend('Mode1', 'Mode2', 'Mode3', 'Mode4', 'Mode5')
hold off

figure (4)
hold on
plot(Zpos30,InterpModes30)
title('Normal Mode Shapes')
xlabel('Distance (m)')
ylabel('InterpModes30')
legend('Mode1', 'Mode2', 'Mode3', 'Mode4', 'Mode5')
hold off


%% Problem 3


MatsSets(1).E = 27e9;
MatsSets(1).A = 0.0003182;
MatsSets(1).rho = 0;

PD.N = 9;
PD.NodePos = [0,     0,   0;
              2,   1.5,   0;
              2,     3,   0;
              0,   4.5,   0;
             -2,   4.5,   0;
             -4,     3,   0;
             -2,     3,   0;
              0,     3,   0;
              0,   1.5,   0];

PD.NE = 15;
PD.ElmConnect = [1, 2;
                 1, 9;
                 2, 3;
                 2, 9;
                 3, 4;
                 3, 8;
                 3, 9;
                 4, 5;
                 4, 8;
                 5, 6;
                 5, 7;
                 5, 8;
                 6, 7;
                 7, 8;
                 8, 9];
PD.NM = 1;
PD.MatsSets = MatsSets;
PD.ElmMats = [1;
              1;
              1;
              1;
              1;
              1;
              1;
              1;
              1;
              1;
              1;
              1;
              1;
              1;
              1];
PD.BCType = [1, 1, 1;
            0, 0, 1;
            0, 0, 1;
            0, 0, 1;
            0, 0, 1;
            0, 1, 1;
            0, 0, 1;
            0, 0, 1;
            0, 0, 1];
PD.BCVal = [0, 0, 0;
       -20000, 0, 0;
       -75000, 0, 0;
            0, -50000, 0;
            0, -40000, 0;
            0, 0, 0;
            0, 0, 0;
            0, 0, 0;
            0, 0, 0];

PDans = PD_truss_static(PD)

NodalDisplacement = PDans.U
MemberStress = PDans.ElmStress

%Truss Displacement
PlotTruss(PDans,1,'y','y')

%% Problem 5

MatsSets(1).E = 210e9;
MatsSets(1).A = 9e-4;
MatsSets(1).rho = 0;

PD2.N = 7;
PD2.NodePos = [0,     0,   0;
               0,   1.5,   0;
            -2.5,   1.5,   0;
             -.5, 0.375,   3;
             -.5, 1.125,   3;
           -1.75, 1.125,   3;
              -1,  0.75,   6];


PD2.NE = 15;
PD2.ElmConnect = [1, 2;
                 1, 3;
                 1, 4;
                 2, 3;
                 2, 4;
                 2, 5;
                 3, 4;
                 3, 5;
                 3, 6;
                 4, 5;
                 4, 6;
                 4, 7;
                 5, 6;
                 5, 7;
                 6, 7];
PD2.NM = 1;
PD2.MatsSets = MatsSets;
PD2.ElmMats = [1;
              1;
              1;
              1;
              1;
              1;
              1;
              1;
              1;
              1;
              1;
              1;
              1;
              1;
              1];
PD2.BCType = [1, 0, 1;
            0, 0, 1;
            1, 1, 1;
            0, 0, 0;
            0, 0, 0;
            0, 0, 0;
            0, 0, 0];

PD2.BCVal = [0, 0, 0;
             0, 0, 0;
             0, 0, 0;
             0, 0, 0;
             0, 0, 0;
             0, 0, 0;
             0, 4000, 0];


PDans2 = PD_truss_static(PD2)

NodalDisplacement2 = PDans2.U
MemberStress2 = PDans2.ElmStress
MaximumMemberStress_On_AD_and_DG = 1.8159
MinimumMemberStress_On_BE_and_EG = -1.8159


%Truss Displacement
PlotTruss(PDans2,1,'y','y')
view([-35.7,30])
