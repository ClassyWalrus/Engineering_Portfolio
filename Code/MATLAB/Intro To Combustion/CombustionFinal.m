%% Combustion Final Project
% Colton Davies
close all; clear; clc;
set(groot,'defaultLineLineWidth',2.0)
set(groot,'defaultAxesXGrid','on')
set(groot,'defaultAxesYGrid','on')
% Well Stirred Reactor Approximation
%Given Parameters
P0 = 101325; %initial pressure - Pa
T0 = 273.15; %K
phi = 0.17;
Re_max = 2000; % Maximum Reynolds number for laminar flow
rhoAir = 1.225; %kg/m^3
rhoCH4 = 0.657; %kg/m^3
rhoC2H6 = 1.35; %kg/m^3
muCH4 = 10.24e-6; %Ns/m^2
muC2H6 = 8.62e-6; %Ns/m^2
muAir = 1.729e-5; %kg/m*s
QFuel = 50*0.00047194745;  % Fuel flow rate in m^3/s
ppm = 1e6;
NOx_limit = 5; % 5 ppm in dry basis

%Calculated Parameters
mdotFuel = QFuel*0.5*rhoCH4+QFuel*0.5*rhoC2H6; %kg/s
QAir = QFuel*5.5/(phi*0.209); %m^3/s
mdotAir = QAir*1.293; %kg/s
mdot = mdotFuel+mdotAir; % kg/s
rhoMix = ((QFuel*0.5*rhoCH4+QFuel*0.5*rhoC2H6)+(QAir*rhoAir))/(QAir+QFuel);
muMix = ((QFuel*0.5*muCH4+QFuel*0.5*muC2H6)+(QAir*muAir))/(QAir+QFuel);

%Nozzle Calcultion
nozzleD = Re_max*muMix/(QFuel*rhoMix); %m
nozzleA = 3.14*nozzleD^2/4;
FuelV = QFuel/nozzleA; %m/s
Re_new = rhoMix*nozzleD*(FuelV)/muMix;
NozD_Final = Re_new/rhoMix/FuelV*muMix;

fprintf('Calculated Nozzle Diameter: %.4f meters\n', nozzleD);
fprintf('Required Combustion Air Flow Rate: %.2f kg/s\n', mdotAir);

destruction_efficiency = 0.99;
residenceTime = 22; %s
reactorVolume = (mdot/rhoMix)*residenceTime
reactorPressure = P0;
InnerD = 2*nozzleD; %m
SteelH = reactorVolume/pi/(InnerD^2);
WallT = 0.003175; %1/8inches in m
OuterD = InnerD+(2*WallT);
TopCap = pi*(InnerD^2)*WallT - pi*(nozzleD^2)*WallT;
SteelVol = pi*SteelH*(OuterD^2-InnerD^2)/4 + TopCap
rho_Steel = 7930; %kg/m^3
SteelMass = rho_Steel*SteelVol %kg
SteelCost = 7.7162*SteelMass % $/kg

%Initial Each Gas Solution
gasFeed = Solution('gri30.yaml');
gasReac = Solution('gri30.yaml');
gasExh = Solution('gri30.yaml'); %this creates a gas phase solution object with 
% species and thermdynamic properties found in the GRI 3.0 mechanism for 
% combustion of ethane and methane.

set(gasFeed,'Temperature',T0,'Pressure',P0);%set initial temperature and pressure of gas

iCH4 = speciesIndex(gasFeed,'CH4'); %indices of species of interest
iC2H6 = speciesIndex(gasFeed, 'C2H6');
iO2  = speciesIndex(gasFeed,'O2');
iN2  = speciesIndex(gasFeed,'N2');
iNO = speciesIndex(gasFeed,'NO');
iOH = speciesIndex(gasFeed, 'OH');

x = zeros(nSpecies(gasFeed),1);
x(iCH4,1) = phi; %set equivalence ratio of initial gas mixture
x(iC2H6,1) = phi;
x(iO2,1) = 5.5;
x(iN2,1) = 3.76*5.5;
set(gasFeed,'Temperature',T0,'Pressure',P0)
set(gasFeed,'MoleFractions',x)
xInit = x;


% x = zeros(nSpecies(gas),1);
x(iCH4,1) = phi; %set equivalence ratio of initial gas mixture
x(iC2H6,1) = phi;
x(iO2,1) = 5.5;
x(iN2,1) = 3.76*5.5;
x(iOH,1) = 5;
set(gasReac,'Temperature',T0,'Pressure',P0)
set(gasReac,'MoleFractions',x)

% create a reactor, and insert the gas
r = IdealGasReactor(gasReac);
setInitialVolume(r,reactorVolume);

%Define Upstream Reservoir
fuelAirMixtureTank = Reservoir(gasFeed);
Exhaust = Reservoir(gasExh);
pReg = Valve(r,Exhaust);
pressureValveCoefficient = 0.01;
setValveCoeff(pReg,pressureValveCoefficient);
% Define mass flow controller for inlet
mfc = MassFlowController(mdot);
install(mfc,r,fuelAirMixtureTank)
% create a reactor network and insert the reactor:
network = ReactorNet({r});

nSteps = 2000;
tim(nSteps) = 0;
temp(nSteps) = 0;
x(nSteps,8) = 0;
t = 0.0;
dt = 1.0e-6;

for n = 1:nSteps
  t = t + dt;
  advance(network, t);
  tim(n) = time(network);
  temp(n) = temperature(r);
  x(n,1:9) = moleFraction(gasReac,{'CH4','C2H6','N','NO', 'NO2', 'OH', 'O2', 'N', 'CO2'});
  X_ex(n,1:6) = moleFraction(gasExh,{'CH4','C2H6','N','NO', 'NO2', 'OH'});
end

%Max concentrations to check reactor function
ReactorMaxT = max(temp(n));
CH4Concentration = max(x(:,1));
C2H6Concentration = max(x(:,2));
NOConcentration = max(x(:,4));
NO2Concentration = max(x(:,5));
OHConcentration = max(x(:,6));
O2Concentration = max(x(:,7));
NConcentration = max(x(:,8));
O2Percent = x(:,7)./(x(:,1)+x(:,2)+...
    x(:,8)+x(:,4)+x(:,6)+x(:,5));
O2Percent = x(:,7)*100-1.5;

%NOx Calculation
NOxConcentration = x(:,4).*ppm+x(:,5).*ppm;
NOx_15O2 = ((NOxConcentration).*(20.9-15))./(20.9-O2Percent);

if max(NOx_15O2)<NOx_limit
    fprintf('NOx Concentration Satisfied: %.2e ppm\n', max(NOx_15O2));
else
    fprintf('NOx Concentration Not Satisfied: %.2e ppm\n', max(NOx_15O2));  
end

%Fuel values for destruction efficiency
CH4 = x(:,1);
CH4max = CH4Concentration;
C2H6 = x(:,2);
C2H6max = C2H6Concentration;

%Destruction Efficiency
D_eff = 1 - (CH4 + C2H6)/(CH4max+C2H6max);

figure(1)
%Temperature
subplot(3,3,1);
plot(tim,temp);
xlabel('Time (s)');
ylabel('Temperature (K)');
%Concentrations
subplot(3,3,2)
semilogx(tim,x(:,1).*ppm);
xlabel('Time (s)');
ylabel('CH4 (ppm)');
subplot(3,3,3)
semilogx(tim,x(:,2).*ppm);
xlabel('Time (s)');
ylabel('C2H6 (ppm)');
subplot(3,3,4)
plot(tim,x(:,4).*ppm);
xlabel('Time (s)');
ylabel('NO (ppm)');
subplot(3,3,5)
plot(tim,x(:,5).*ppm);
xlabel('Time (s)');
ylabel('NO2 (ppm)');
subplot(3,3,6)
semilogx(tim,x(:,9).*ppm);
xlabel('Time (s)');
ylabel('CO2 (ppm)');
%Destruction
subplot(3,3,7)
semilogx(tim,D_eff);
xlabel('Time (s)');
ylabel('Destruction Efficiency %');
%NOx
subplot(3,3,8)
semilogx(tim,NOxConcentration);
hold on
yline(5, 'r')
xlabel('Time (s)');
ylabel('NOx (ppm)');
subplot(3,3,9)
semilogx(tim,NOx_15O2);
hold on
yline(5, 'r')
xlabel('Time (s)');
ylabel('NOx at 15% O2 (ppm)');

figure(3)
semilogx(tim,x(:,7).*ppm);
xlabel('Time (s)');
ylabel('O2 (ppm)');


%% Premixed Free Flame

% FLAME1 - A burner-stabilized flat flame
%
%    This script simulates a burner-stablized lean hydrogen-oxygen flame
%    at low pressure.
%
% Keywords: combustion, 1D flow, burner-stabilized flame, plotting

% help flame1

t0 = cputime;  % record the starting time

% parameter values
p          =   0.05*oneatm;         % pressure
tburner    =   273.15;               % burner temperature

rxnmech    =  'gri30.yaml';           % reaction mechanism file
comp       =  'CH4:0.17, C2H6:0.17, N2:20.68, O2:5.5, OH:5'; % premixed gas composition
mdot = mdot/nozzleA;
initial_grid = [0.0 0.02 0.04 0.06 0.08 0.1 ...
                0.15 0.2 0.4 0.49 0.5];  % m

initial_grid = linspace(0,SteelH,100);

tol_ss    = [1.0e-4 1.0e-10];       % [rtol atol] for steady-state
                                    % problem
tol_ts    = [1.0e-3 1.0e-7];        % [rtol atol] for time stepping

loglevel  = 0;                      % amount of diagnostic output (0
                                    % to 5)
refine_grid = 1;                    % 1 to enable refinement, 0 to
                                    % disable
max_jacobian_age = [5, 10];

%%%%%%%%%%%%%%%% create the gas object %%%%%%%%%%%%%%%%%%%%%%%%
%
% This object will be used to evaluate all thermodynamic, kinetic,
% and transport properties
%
% gas = Solution(rxnmech, 'ohmech', 'mixture-averaged');

% set its state to that of the unburned gas at the burner
set(gasFeed,'T', tburner, 'P', p, 'X', comp);

%%%%%%%%%%%%%%%% create the flow object %%%%%%%%%%%%%%%%%%%%%%%

f = AxisymmetricFlow(gasFeed,'flow');
set(f, 'P', p, 'grid', initial_grid);
set(f, 'tol', tol_ss, 'tol-time', tol_ts);

%%%%%%%%%%%%%%% create the burner %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%  The burner is an Inlet object. The temperature, mass flux,
%  and composition (relative molar) may be specified.
%
burner = Inlet('burner');
set(burner, 'T', tburner, 'MassFlux', mdot, 'X', comp);

%%%%%%%%%%%%%% create the outlet %%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%  
%
s = Outlet('out');

%%%%%%%%%%%%% create the flame object  %%%%%%%%%%%%
%
% Once the component parts have been created, they can be assembled
% to create the flame object.
%
fl = flame(gasFeed, burner, f, s);
setMaxJacAge(fl, max_jacobian_age(1),  max_jacobian_age(2));

% if the starting solution is to be read from a previously-saved
% solution, uncomment this line and edit the file name and solution id.
%restore(fl,'h2flame.yaml', 'energy')

solve(fl, loglevel, refine_grid);

%%%%%%%%%%%% enable the energy equation %%%%%%%%%%%%%%%%%%%%%
%
%  The energy equation will now be solved to compute the
%  temperature profile. We also tighten the grid refinement
%  criteria to get an accurate final solution.
%

enableEnergy(f);
setRefineCriteria(fl, 1, 200.0, 0.1, 0.1);
solve(fl, 1, 0.1);
% saveSoln(fl,'gri30.yaml','energy',['solution with energy equation']);

%%%%%%%%%% show statistics %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
writeStats(fl);
elapsed = cputime - t0;
e = sprintf('Elapsed CPU time: %10.4g',elapsed);
disp(e);

%%%%%%%%%% make plots %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure(2)
subplot(3,3,1);
plotSolution(fl, 'flow', 'T');
title('Temperature [K]');
subplot(3,3,2);
plotSolution(fl, 'flow', 'velocity');
title('Axial Velocity [m/s]');
subplot(3,3,3);
plotSolution(fl, 'flow', 'O2');
title('O2 Mass Fraction');
subplot(3,3,4);
plotSolution(fl, 'flow', 'NO');
title('NO Mass Fraction');
subplot(3,3,5);
plotSolution(fl, 'flow', 'NO2');
title('NO2 Mass Fraction');
subplot(3,3,6);
plotSolution(fl, 'flow', 'H2O');
title('H2O Mass Fraction');
subplot(3,3,7);
plotSolution(fl, 'flow', 'CH4');
title('CH4 Mass Fraction');
subplot(3,3,8);
plotSolution(fl, 'flow', 'C2H6');
title('C2H6 Mass Fraction');
subplot(3,3,8);
plotSolution(fl, 'flow', 'CO2');
title('CO2 Mass Fraction');


%% Laminar Diffusion using Roper Flame Height Correlation
S = 17.19;
Tinf = 273.15; %K
Tfuel = 1500; %K
Lf = 1330*(QFuel*(Tinf/Tfuel))/(log(1+1/S)) %m
Ratio_Laminar = Lf/nozzleD


%% Turbulent Non-premixed flame
Re_Turbulent = 20000;
Nozzle_Turbulent = Re_Turbulent*muMix/(QFuel*rhoMix) %m
Ratio_Turbulent = 108
Lf_TurbulentEstimate = Ratio_Turbulent*Nozzle_Turbulent


%% CO2 Emmissions
VolumetricFlowRate = 0.02359735999901; %m^3/s
Seconds = 3.154e7; %s
VolumePerYear = VolumetricFlowRate*Seconds;
CH4_GWP = 30;
C2H6_GWP = 10;
mdot_CH4 = VolumetricFlowRate/2*rhoCH4*Seconds; %kg/year
mdot_C2H6 = VolumetricFlowRate/2*rhoC2H6*Seconds; %kg/year

CO2_EmissionsCH4 = mdot_CH4*2.75 %kg/year
CO2_EmissionsC2H6 = mdot_C2H6*2.93 %kg/year

UnburntCH4_Equivalent = 30*CO2_EmissionsCH4 %kg/year
UnburntC2H6_Equivalent = 10*CO2_EmissionsC2H6 %kg/year

%From molar mass ratio, we know that:
    %1kg of methane produces 2.75kg of CO2
    %1kg of ethane produces 2.93kg of CO2

