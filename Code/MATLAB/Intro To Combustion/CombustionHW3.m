clc;
clear;

g = GRI30;
% set(gas1,'T',900.0,'P',1.e5,'X','CH4:1,O2:2,N2:7.52');
disp('fixed U and V:');
set(g,'T',1200.0,'P',oneatm,'X','CH4:0.95,O2:2,N2:7.52');
equilibrate(g,'UV')
disp('fixed S and V:');
set(g,'T',1200.0,'P',oneatm,'X','CH4:0.95,O2:2,N2:7.52');
equilibrate(g,'SV')
disp('fixed S and P:');
set(g,'T',1200.0,'P',oneatm,'X','CH4:0.95,O2:2,N2:7.52');
equilibrate(g,'SP')

