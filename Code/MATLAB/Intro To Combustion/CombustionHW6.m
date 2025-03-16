%% Combustion HW6
% Colton Davies

%% Problem 2

 A = 1.06036; B = 0.15610;
 C = 0.19300; D = 0.47635;
 E = 1.03587; F = 1.52996;
 G = 1.76474; H = 3.89411;

k_B = 1.380649e-23;
MW_He= 4.0026;
MW_O2 = 32;
MW_CH4 = 16.04;

MW_HeO2 = 2*(1/MW_He+1/MW_O2)^(1/2);
MW_HeCH4 = 2*(1/MW_He+1/MW_CH4)^(1/2);
MW_O2CH4 = 2*(1/MW_O2+1/MW_CH4)^(1/2);

T = 500;
P = 101325;

%Data from table D
sig_He = 2.551;  
sig_O2 = 3.467;
sig_CH4 = 3.758;
ep_He = 2.551*k_B;
ep_O2 = 106.7*k_B;
ep_CH4 = 148.6*k_B;


sig_HeO2 = (sig_He+sig_O2)/2;
sig_HeCH4 = (sig_He+sig_CH4)/2;
sig_O2CH4 = (sig_O2+sig_CH4)/2;

T_star_HeO2 = k_B*T/(ep_He*ep_O2)^1/2;
T_star_HeCH4 = k_B*T/(ep_He*ep_CH4)^1/2;
T_star_O2CH4 = k_B*T/(ep_CH4*ep_O2)^1/2;

Omega_D_HeO2 = A/T_star_HeO2^B+C/exp(D*T_star_HeO2)+E/exp(F*T_star_HeO2)+G/exp(H*T_star_HeO2);
Omega_D_HeCH4 = A/T_star_HeCH4^B+C/exp(D*T_star_HeCH4)+E/exp(F*T_star_HeCH4)+G/exp(H*T_star_HeCH4);
Omega_D_O2CH4 = A/T_star_O2CH4^B+C/exp(D*T_star_O2CH4)+E/exp(F*T_star_O2CH4)+G/exp(H*T_star_O2CH4);

% Solve for individual effective binary diffusion coefficients
D_HeO2 = (0.0266*T^(3/2))/(P*MW_HeO2^(1/2)*sig_HeO2^(2)*Omega_D_HeO2)
D_HeCH4 = (0.0266*T^(3/2))/(P*MW_HeCH4^(1/2)*sig_HeCH4^(2)*Omega_D_HeCH4)
D_O2CH4 = (0.0266*T^(3/2))/(P*MW_O2CH4^(1/2)*sig_O2CH4^(2)*Omega_D_O2CH4)

Deff_He = (1 - 1/3)/((1/3)/D_HeO2+(1/3)/D_HeCH4)
Deff_O2 = (1 - 1/3)/((1/3)/D_HeO2+(1/3)/D_O2CH4)
Deff_CH4 = (1 - 1/3)/((1/3)/D_HeCH4+(1/3)/D_O2CH4)


%% Problem 3

gas = Solution('gri30.yaml');
P0 = 1e5;%initial pressure (Pa)
T0 = 1000; %initial temperature (K)
nFs = 50;
F = linspace(0,0.06,nFs);
tad = zeros(1,nFs);

%State Relation for mixture fraction to phi
phi = (F./(1-F))/(0.0527);


xResults = zeros(nSpecies(gas),nFs); 
% set_mixture_fraction(0.5, 'CH4', 'O2:1.0, N2:3.76')
% mixture_fraction(fuel, oxidizer, basis='mole', element='Bilger')
for i = 1:1:nFs
    gas = Solution('gri30.yaml'); %this creates a gas phase solution object with 
    % species and thermdynamic properties found in the GRI 3.0 mechanism for 
    % combustion of methane.

    set(gas,'Temperature',T0,'Pressure',P0);%set initial temperature and pressure of gas
% set(mixture_fraction(F(i), 'CH4', 'O2:1.0, N2:3.76'))
    iCH4 = speciesIndex(gas,'CH4'); %indices of species of interest
    iO2  = speciesIndex(gas,'O2');
    iN2  = speciesIndex(gas,'N2');
    iNO = speciesIndex(gas,'NO');

    x = zeros(nSpecies(gas),1);
    x(iCH4,1) = phi(i); %set equivalence ratio of initial gas mixture
    x(iO2,1) = 2;
    x(iN2,1) = 3.76*2;

    set(gas,'MoleFractions',x)
    equilibrate(gas,'HP'); %calculates equilibrium conditions under fixed specific enthalpy and pressure

    tad(i) = temperature(gas); %equilibrium temperature (K)
    xResults(:,i) = moleFractions(gas); %species mole fractions at equilibrium

end

figure(1)
hold on
plot(F,tad,marker="x");
title('Adiabatic Flame Temp for Mixture Fraction')
xlabel('Mixture Fraction');
ylabel('Max Temperature (K)');
legend('Equilibrate')
grid on



