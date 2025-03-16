%%Combustion HW10
%Colton Davies
close all; clear; clc;

Ts = [500, 1000, 1500, 2000];
MW_C = 0.01201;  %kg/kmol
p =  2e17;  %kg/m^3

%For 1kg of C and V=1m^3
Conc_C_1 = 0.083257;  %kmol/m^3 0.083257  1341.5566
Conc_C_2 = 0.083257;

k1 = 3.007e7*exp(-17966./Ts);  %m/s
k2 = 4.016e8*exp(-29790./Ts);  %m/s

R1 = k1*(Conc_C_1^2)*MW_C;  %kg_C/m^2*s
R2 = k2*(Conc_C_2)*MW_C;  %kg_C/m^2*s

Ratio = R1./R2;

R1_Values = array2table(R1, 'VariableNames', {'500K', '1000K', '1500K', '2000K'})
R2_Values = array2table(R2, 'VariableNames', {'500K', '1000K', '1500K', '2000K'})
Reaction_Rate_Ratio = array2table(Ratio, 'VariableNames', {'500K', '1000K', '1500K', '2000K'})