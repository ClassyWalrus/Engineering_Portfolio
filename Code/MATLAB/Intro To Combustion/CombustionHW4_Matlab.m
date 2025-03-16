% clc;clear;close all;
% %Example code for in-class work on 09-17-2024
% 
% %Author: Cliff Goertemiller
% %Date of Creation: 09-15-2024
% 
% %The purpose of this code is to show students how to perform the following
% %tasks in Cantera
% %1. create a gas phase solution object with properties based on a chemical
% %kinetic mechanism
% %2. set the initial pressure and temperature of that gas phase solution
% %object
% %3. set up data structures to store the equivalence ratios of interest the 
% % and the equilibrium species concentration results
% %4. calculate an equilibrium concentration for a combustable gas mixture at
% % each equivalance ratio based on constant specific enthalpy and constant
% % volume (no heat or work leaves the system)
% %5. calculate the adiabatic flame temperature based on chemical equilibrium
% %of each gas mixture under constant specific enthalpy and constant pressure
% %conditions
% %6. present results on a figure
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% gas = Solution('gri30.yaml');
% P0 = 1e5;%initial pressure (Pa)
% T0 = 1000; %initial temperature (K)
% P=oneatm;
% 
% nPhis = 25; %number of equivalence ratios to evaluate
% 
% phi = linspace(0.6,1.6,nPhis);
% tad = zeros(1,nPhis);
% xResults = zeros(nSpecies(gas),nPhis); 
% 
% for i = 1:1:nPhis
%     gas = Solution('gri30.yaml'); %this creates a gas phase solution object with 
%     % species and thermdynamic properties found in the GRI 3.0 mechanism for 
%     % combustion of methane.
% 
%     set(gas,'Temperature',T0,'Pressure',P0);%set initial temperature and pressure of gas
% 
%     iCH4 = speciesIndex(gas,'CH4'); %indices of species of interest
%     iO2  = speciesIndex(gas,'O2');
%     iN2  = speciesIndex(gas,'N2');
%     iNO = speciesIndex(gas,'NO');
% 
%     x = zeros(nSpecies(gas),1);
%     x(iCH4,1) = phi(i); %set equivalence ratio of initial gas mixture
%     x(iO2,1) = 2;
%     x(iN2,1) = 3.76*2;
% 
%     set(gas,'MoleFractions',x)
%     equilibrate(gas,'HP'); %calculates equilibrium conditions under fixed specific enthalpy and pressure
% 
%     tad(i) = temperature(gas); %equilibrium temperature (K)
%     xResults(:,i) = moleFractions(gas); %species mole fractions at equilibrium
% 
% 
% % iNO = speciesIndex(gas,'NO');
% % iNO2 = speciesIndex(gas,'NO2');
% % iCO = speciesIndex(gas,'CO');
% % iH2O = speciesIndex(gas,'H2O');
% % iCO2 = speciesIndex(gas,'CO2');
% 
% 
% % create a reactor, and insert the gas
% r = Reactor(gas);
% 
% % create a reservoir to represent the environment
% a = Solution('air.yaml','air','none');
% set(a,'P',P)
% env = Reservoir(a);
% 
% % Define a wall between the reactor and the environment and
% % make it flexible, so that the pressure in the reactor is held
% % at the environment pressure.
% w = Wall;
% install(w,r,env);
% 
% % set expansion parameter. dV/dt = KA(P_1 - P_2)
% setExpansionRateCoeff(w, 1.0e6);
% 
% % set wall area
% setArea(w, 1.0);
% 
% % create a reactor network and insert the reactor:
% network = ReactorNet({r});
% 
% 
% 
% nSteps = 2000;
% tim(nSteps) = 0;
% temp(nSteps) = 0;
% x(nSteps,3) = 0;
% t = 0.0;
% dt = 1.0e-6;
% t0 = cputime;
% 
%     for n = 1:nSteps
%       t = t + dt;
%       advance(network, t);
%       tim(n) = time(network);
%       temp(n) = temperature(r);
%       x(n,1:3) = moleFraction(gas,{'CH4','N','NO'});
%     end
% 
%     ReactorMaxT(i) = max(temp(n));
%     NOConcentration(i) = max(x(:,3));
% 
% end
% 
% 
% 
% % figure(1)
% % hold on;
% % plot(phi,tad)
% % xlabel('\phi')
% % ylabel('Equilibrium Temperature (K)')
% % grid on;
% % 
% % figure(2)
% % hold on;
% % plot(phi,xResults(iNO,:).*1e6)
% % % plot(phi,xResults(iCO,:).*1e6)
% % % plot(phi,xResults(iO2,:).*1e6)
% % % plot(phi,xResults(iH2O,:).*1e6)
% % % plot(phi,xResults(iCO2,:).*1e6)
% % xlabel('\phi')
% % ylabel('Equilibrium Cocentration (PPM)')
% % 
% % yyaxis right
% % hold on;
% % % plot(phi,xResults(iNO2,:).*1e6)
% % plot(phi,xResults(iCH4,:).*1e6)
% % 
% % grid on;
% % legend('NO','CH_4')
% % title('Species Concentrations')
% % hold off
% % 
% % figure(3)
% % grid on
% % subplot(2,2,1);
% % plot(tim,temp);
% % xlabel('Time (s)');
% % ylabel('Temperature (K)');
% % subplot(2,2,2)
% % plot(tim,x(:,1));
% % xlabel('Time (s)');
% % ylabel('CH4 Mole Fraction (K)');
% % subplot(2,2,3)
% % plot(tim,x(:,2));
% % xlabel('Time (s)');
% % ylabel('N Mole Fraction (K)');
% % subplot(2,2,4)
% % plot(tim,x(:,3));
% % xlabel('Time (s)');
% % ylabel('NO Mole Fraction (K)');
% % clear all
% % cleanup
% 
% 
% figure(4)
% hold on
% plot(phi,tad,marker="x");
% plot(phi,ReactorMaxT,marker="o");
% title('Equilibrate VS Reactor')
% xlabel('Phi');
% ylabel('Max Temperature (K)');
% legend('Equilibrate','Reactor')
% grid on
% 
% 
% figure(5)
% hold on
% plot(phi,xResults(iNO,:),marker="x");
% plot(phi,NOConcentration,marker="o");
% title('Equilibrate VS Reactor')
% xlabel('Phi');
% ylabel('NO Concentration');
% legend('Equilibrate','Reactor')
% grid on
% 

%% Problem 4
clc; clear; close all;
P0 = 101325; %initial pressure - Pa
T0 = 1000; %K
phi = linspace(0.6,1.6,20);
maxtemp(length(phi)) =0;
ignDelayOut(length(phi)) =0;

for k=1:1:length(phi)
gas = Solution('gri30.yaml');
set(gas,'Temperature',T0,'Pressure',P0);%set initial temperature and pressure of gas
iCH4 = speciesIndex(gas,'CH4'); %indices of species of interest
iO2 = speciesIndex(gas,'O2');
iN2 = speciesIndex(gas,'N2');
x = zeros(nSpecies(gas),1);
x(iCH4,1) = phi(k); %set equivalence ratio of initial gas mixture
x(iO2,1) = 2;
x(iN2,1) = 3.76*2;
set(gas,'MoleFractions',x)
% set the initial conditions
% create a reactor, and insert the gas
r = IdealGasReactor(gas);
% create a reservoir to represent the environment
a = Solution('air.yaml','air','none');
set(a,'P',P0)
env = Reservoir(a);
% Define a wall between the reactor and the environment and
% make it flexible, so that the pressure in the reactor is held
% at the environment pressure.
w = Wall;
install(w,r,env);
% set expansion parameter. dV/dt = KA(P_1 - P_2)
setExpansionRateCoeff(w, 1.0e6);
% set wall area
setArea(w, 1.0);
% create a reactor network and insert the reactor:
network = ReactorNet({r});
nSteps = 1000;% number of time steps in simulation
tim(nSteps) = 0;
temp(nSteps) = 0;
x(nSteps,3) = 0;
HRR(nSteps) = 0;
t = 0.0; %initial simulation time
dt = 1.0e-6; %time step duration (s)
t0 = cputime;
for n = 1:nSteps
t = t + dt;
advance(network, t);
tim(n) = time(network);
temp(n) = temperature(r);
x(n,1:4) = moleFraction(gas,{'OH','H','H2','NO'});
netProd = netProdRates(gas);
HRR(n) = sum(netProd.*enthalpy_mole(gas));
end
maxTemp(k) = max(temp);
[~,idxIgnDelay]=max(x(:,1));
ignDelayOut(k) = tim(idxIgnDelay);
maxHRR(k) = max(HRR);
maxNO(k) = max(x(:,4));
end
%disp(['CPU time = ' num2str(cputime - t0)]);
clf;
subplot(2,2,1);
plot(tim,temp);
xlabel('Time (s)');
ylabel('Temperature (K)');
subplot(2,2,2)
plot(tim,x(:,1));
xlabel('Time (s)');
ylabel('OH Mole Fraction (K)');
subplot(2,2,3)
plot(tim,x(:,2));
xlabel('Time (s)');
ylabel('H Mole Fraction (K)');
subplot(2,2,4)
plot(tim,x(:,3));
xlabel('Time (s)');
ylabel('H2 Mole Fraction (K)');
%now calculate the ignition delay
%igniton delay is the duration between the start of combustion and the time
%at which OH concentrations peak
disp("The Igniton delay is " + tim(idxIgnDelay) + " seconds")
figure(2)
plot(phi,maxHRR)
xlabel('Time (s)')
ylabel('maximum HRR (W/m^3)')
figure(3)
plot(phi,maxTemp)
xlabel('\phi')
ylabel('Maximum Temperature (K)')
figure(4)
plot(phi,ignDelayOut.*1000)
xlabel('\phi')
ylabel('Ignition Delay (ms)')
figure(5)
plot(phi,maxNO*1e6)
xlabel('\phi')
ylabel('peak NO concentration (PPM)')
cleanup