% MATLAB code to simulate a 1D flat H2-O2 flame diluted with argon

% Load Cantera
% addpath('path_to_cantera_matalb_api'); % Make sure to add Cantera's Matlab interface

% Create the gas mixture using the 'gri30' mechanism (which includes H2, O2, Ar)
gas = Solution('h2o2.yaml');

% Define the molar ratios for H2, O2, and Ar
H2 = 2.0;  % 2 moles of H2 per mole of fuel
O2 = 1.0;  % 1 mole of O2
Ar = 7.0;  % 7 moles of Ar per mole of O2

% Equivalence ratio
phi = 0.9;

% Calculate the fuel-to-oxidizer ratio for the given equivalence ratio
% Using phi = actual fuel/oxidizer / stoichiometric fuel/oxidizer ratio
stoich_fuel_O2_ratio = 2.0; % Stoichiometric H2:O2 ratio is 2:1
fuel = phi * O2 * stoich_fuel_O2_ratio / (1 + phi * stoich_fuel_O2_ratio);
oxidizer = O2 / (1 + phi * stoich_fuel_O2_ratio);
argon = Ar * O2 / (O2 + Ar);

% Set the gas mixture composition based on the equivalence ratio
setEquivalenceRatio(phi, fuel, oxidizer, argon);

% Set the initial conditions (temperature, pressure, and gas composition)
TPX = [1200, 101325, gas.moleFractions];  % Initial temperature 1200 K, 1 atm pressure

% Define the 1D flame problem
% Create a FreeFlame object
flame = FreeFlame(gas);

% Set the grid size (this will control the resolution of the solution)
setGrid(0.001, 0.1);  % Flame width between 0.001 and 0.1 meters

% Solve the flame structure
solve(flame);

% Extract the results
x = flame.grid;  % Distance along the flame front
T = flame.T;     % Temperature profile (K)
O2 = flame.Y(:, gas.speciesIndex('O2'));  % O2 mole fraction profile

% Plot the temperature profile
figure;
subplot(2,1,1);
plot(x, T, 'LineWidth', 2);
xlabel('Distance (m)');
ylabel('Temperature (K)');
title('Temperature Profile');
grid on;

% Plot the O2 mole fraction profile
subplot(2,1,2);
plot(x, O2, 'LineWidth', 2, 'Color', 'orange');
xlabel('Distance (m)');
ylabel('O2 Mole Fraction');
title('O2 Profile');
grid on;

% Adjust the layout of the plots
tight_layout;