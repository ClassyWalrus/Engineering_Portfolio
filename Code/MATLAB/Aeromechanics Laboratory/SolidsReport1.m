%% Aeromechanics Solids Report
% Colton Davies
close all; clear; clc;
set(groot,'defaultLineLineWidth',1.5)

%% Extensometer Calibration
%Channel 1 = Force % Tensile Machine:562lbs/Volt || DIC Tensile:690lbs/Volt
%Channel 2 = Crosshead position % 2in = 10V
%Channel 3 = Extensometer % 2in = 10V
    in2mm = 25.4; %inches to mm
    lbs2N = 4.44822; %pounds to newtons
    Calibration_Pos_Extensometer = [-.1 -0.05 0.05 0.1 0.15 0.2]*in2mm;
    Calibration_V = [-0.95 -0.474 0.47 0.948 1.424 1.898];
    Extensometer_LinearRegression = polyfit(Calibration_V, Calibration_Pos_Extensometer, 1);
    %Instrumentation Conversion Factors
    Extensometer_Conv = Extensometer_LinearRegression(1); %%mm/volt
    Crosshead_Conv = 50.8/10; %mm/volt
    DIC_Conv = 690*lbs2N; %N/volt
    TensileTest_Conv = 562*lbs2N; %N/volt
    %Extensometer removed indices
    IDX_S_E_Removal = 434-62;
    IDX_A_E_Removal = 222-70;
    IDX_DIC_E_Removal = 140;
    %Fracture indices
    IDX_S_Fracture = 578-62;
    IDX_A_Fracture = 230-70;
    IDX_DIC_Fracture = 145;
%% Material Properties and Measurements
%Stainless Steel
    %Shape parameters
    t_S_init = 0.056*in2mm;
    t_S_fin = t_S_init-.1;
    w_S_init = 0.405*in2mm;
    w_S_fin = w_S_init-3.3;
    GaugeL_S = 2*in2mm;
    A0_S = w_S_init*t_S_init; % Initial Cross Section area mm^2
    %import csv
    S_Data = readmatrix('ss_1008_test.csv');
    S_time = S_Data(2:end,1);
    %Force and Extension
    S_Force = S_Data(62:end,2)*TensileTest_Conv; % Newton Force from Volt Conversion
    S_Crosshead = S_Data(62:end,4)*Crosshead_Conv;
    S_Extensometer = S_Data(62:end,6)*Extensometer_Conv;

%Aluminum 6061
    %Shape parameters
    t_A_init = 0.062*in2mm;
    t_A_fin = t_A_init-0.1;
    w_A_init = 0.405*in2mm;
    w_A_fin = w_A_init-0.4;
    GaugeL_A = 2*in2mm;
    A0_A = w_A_init*t_A_init; %Initial Cross Section area mm^2
    %import csv
    A_Data = readmatrix('ai_6061.csv');
    A_time = A_Data(61:end,1);
    %Force and Extension
    A_Force = A_Data(70:end,2)*TensileTest_Conv; % Newton Force from Volt Conversion
    A_Crosshead = A_Data(70:end,4)*Crosshead_Conv; % mm of extension
    A_Extensometer = A_Data(70:end,6)*Extensometer_Conv; % mm of extension

%DIC
    % Subset size 33, step size 8, ssig threshold 150, gauge size 24
    % Images used(32,74,102,123,143)
    %Shape parameters
    t_DIC_init = 0.062*in2mm; %mm
    t_DIC_fin = t_DIC_init-0.1; %mm
    w_DIC_init = 0.601*in2mm; %mm
    w_DIC_fin = w_DIC_init-2; %mm
    w_DIC_hole_init = w_DIC_init/3; %mm
    w_DIC_hole_fin = w_DIC_hole_init / (69/81); %mm
    GaugeL_DIC = 2*in2mm; %mm
    A0_DIC = (w_DIC_init-w_DIC_hole_init)*t_DIC_init; %Inital Cross Sectional area mm^2
    %import csv
    DIC_Data = readmatrix('ss_hole.csv');
    DIC_Data(1,2) = 0;
        DIC_Data(1,4) = 0;
            DIC_Data(1,6) = 0;
    DIC_time = DIC_Data(2:end,1);
    %import DIC analysis data
    DIC_Anyalysis = readmatrix('live_plot_line_frame_0.txt','Delimiter',{','});
    VSG_STRAIN_XX = DIC_Anyalysis(2:end,8);
    VSG_STRAIN_XY = DIC_Anyalysis(2:end,9);
    VSG_STRAIN_YY = DIC_Anyalysis(2:end,10);
    ARC_LENGTH = DIC_Anyalysis(2:end,1); %pixels
    %Force and Extension
    DIC_Force = DIC_Data(1:end,2)*DIC_Conv; % Force from Volt Conversion
    DIC_Crosshead = DIC_Data(1:end,4)*Crosshead_Conv;
    DIC_Extensometer = DIC_Data(1:end,6)*Extensometer_Conv;
    
%% Combine Crosshead and Extensometer Data
%Steel
    L_S = zeros(IDX_S_Fracture,1);
    L_S(1:IDX_S_E_Removal) = GaugeL_S + S_Extensometer(1:IDX_S_E_Removal);
    L_S(IDX_S_E_Removal-1:IDX_S_Fracture) = ...
        L_S(IDX_S_E_Removal-1) + ... % Last extensometer point
        (S_Crosshead(IDX_S_E_Removal-1:IDX_S_Fracture) - ...
        S_Crosshead(IDX_S_E_Removal-1)); % Crosshead displacement
%Aluminum
    L_A = zeros(IDX_A_Fracture,1);
    L_A(1:IDX_A_E_Removal) = GaugeL_A + ...
        A_Extensometer(1:IDX_A_E_Removal);
    L_A(IDX_A_E_Removal-1:IDX_A_Fracture) = ...
        L_A(IDX_A_E_Removal-1) + ... % Last extensometer point
        (A_Crosshead(IDX_A_E_Removal-1:IDX_A_Fracture) - ...
        A_Crosshead(IDX_A_E_Removal-1)); % Crosshead displacement
%DIC
    L_DIC = zeros(IDX_DIC_Fracture,1);
    L_DIC(1:IDX_DIC_E_Removal) = GaugeL_DIC + ...
        DIC_Extensometer(1:IDX_DIC_E_Removal);
    L_DIC(IDX_DIC_E_Removal-1:IDX_DIC_Fracture) = ...
        L_DIC(IDX_DIC_E_Removal-1) + ... % Last extensometer point
        (DIC_Crosshead(IDX_DIC_E_Removal-1:IDX_DIC_Fracture) - ...
        DIC_Crosshead(IDX_DIC_E_Removal-1)); % Crosshead displacement

%% Stress and Strain [MPa, Unitless]
%Stainless Steel
    %Engineering Stress and Strain
    S_Eng_Stress = S_Force(1:IDX_S_Fracture)/A0_S;
    S_Eng_Strain = (L_S-GaugeL_S)/GaugeL_S;
    %True Stress and Strain
    S_Instant_Area = A0_S*GaugeL_S./L_S;
    S_True_Stress = S_Force(1:IDX_S_Fracture)./S_Instant_Area;
    S_True_Strain = log(L_S/GaugeL_S);
    %Young's Modulus
    S_Polyfit= polyfit(S_Eng_Strain(1:8), S_Eng_Stress(1:8),1);
    S_Fit= fit(S_Eng_Strain(1:8), S_Eng_Stress(1:8),'poly1');
    S_Youngs = S_Polyfit(1);
    %Poisson Ratio
    PR_S = -VSG_STRAIN_XX/VSG_STRAIN_YY;
    %Yield Stress
    S_Offset_Strain = linspace(-0.002,0.001,578)';
   % S_Offset_Stress = polyval(S_fit(1:2)',S_Offset_Strain);
    S_Offset_Stress = polyval(S_Polyfit(1:2)',S_Offset_Strain);
    [IDX_S_Dist, S_Dist_Min] = dsearchn([S_Offset_Strain + 0.002,...
    S_Offset_Stress], [S_Eng_Strain S_Eng_Stress]);
    [S_Dist_Min, IDX_Dist_Min] = min(S_Dist_Min);
    S_Yield_Stress = S_Offset_Stress(IDX_S_Dist(IDX_Dist_Min));
    %Ultimate Strength
    S_Ult_Stress = max(S_Eng_Stress);
    %Percent Reduction of Area at Failure
    Afin_S = w_S_fin*t_S_fin;
    Red_S = 1-Afin_S/A0_S;
    %True Strain at Failure
    TS_S = S_True_Strain(end);
    %Toughness
    Tough_S = cumtrapz(S_Eng_Strain, S_Eng_Stress);
    Tough_S = Tough_S(end);
%Aluminum 6061
    %Engineering Stress and Strain
    A_Eng_Stress = A_Force(1:IDX_A_Fracture)/A0_A;
    A_Eng_Strain = (L_A-GaugeL_A)/GaugeL_A;
    %True Stress and Strain
    A_Instant_Area = A0_A*GaugeL_A./L_A;
    A_True_Stress = A_Force(1:IDX_A_Fracture)./A_Instant_Area;
    A_True_Strain = log(L_A/GaugeL_A);
    %Young's Modulus
    A_Polyfit = polyfit(A_Eng_Strain(25:50), A_Eng_Stress(25:50),1);
    A_Fit= fit(A_Eng_Strain(25:50), A_Eng_Stress(25:50),'poly1');
    A_Youngs = A_Polyfit(1);
    %Yield Stress
    A_Offset_Strain = linspace(-0.002,0.006,230)';
    A_Offset_Stress = polyval(A_Polyfit,A_Offset_Strain);
    [IDX_A_Dist, A_Dist_Min] = dsearchn([A_Offset_Strain + 0.002,...
    A_Offset_Stress], [A_Eng_Strain A_Eng_Stress]);
    [A_Dist_Min, IDX_Dist_Min] = min(A_Dist_Min);
    A_Yield_Stress = A_Offset_Stress(IDX_A_Dist(IDX_Dist_Min));
    %Ultimate Strength
    A_Ult_Stress = max(A_Eng_Stress);
    %Percent Reduction of Area at Failure
    Afin_A = w_A_fin*t_A_fin;
    Red_A = 1-Afin_A/A0_A;
    %True Strain at Failure
    TS_A = A_True_Strain(end);
    %Toughness
    Tough_A = cumtrapz(A_Eng_Strain, A_Eng_Stress);
    Tough_A = Tough_A(end);
%DIC
    %Engineering Stress and Strain
    DIC_Eng_Stress = DIC_Force(1:IDX_DIC_Fracture)/A0_DIC;
    DIC_Eng_Strain = (L_DIC-GaugeL_DIC)/GaugeL_DIC;
    %True Stress and Strain
    DIC_Instant_Area = A0_DIC*GaugeL_DIC./L_DIC;
    DIC_True_Stress = DIC_Force(1:IDX_DIC_Fracture)./DIC_Instant_Area;
    DIC_True_Strain = log(L_DIC/GaugeL_DIC);
    %Young's Modulus
    DIC_Polyfit = polyfit(DIC_Eng_Strain(3:25), DIC_Eng_Stress(3:25),1);
    DIC_Fit= fit(DIC_Eng_Strain(3:25), DIC_Eng_Stress(3:25),'poly1');
    DIC_Youngs = DIC_Polyfit(1);
    %Yield Stress
    DIC_Offset_Strain = linspace(-0.002,0.001,145)';
    DIC_Offset_Stress = polyval(DIC_Polyfit,DIC_Offset_Strain);
    %Ultimate Strength
    DIC_Ult_Stress = max(DIC_Eng_Stress);
    %Percent Reduction of Area at Failure
    Afin_DIC = (w_DIC_fin-w_DIC_hole_fin)*t_DIC_fin;
    Red_DIC = 1-Afin_DIC/A0_DIC;
    %True Strain at Failure
    TS_DIC = DIC_True_Strain(end);
    %Toughness
    Tough_DIC = cumtrapz(DIC_Eng_Strain, DIC_Eng_Stress);
    Tough_DIC = Tough_DIC(end);
%%
%% Plotting
%Plotting Steel
    %Stress and Strain
    figure(1)
    hold on
    plot(S_Eng_Strain, S_Eng_Stress, "Color","k","LineStyle","-")
    plot(S_True_Strain, S_True_Stress, "Color","r","LineStyle","--")
    title("1008 Steel Stress-Strain Curves")
    xlabel("Strain $\epsilon$ [mm/mm]","Interpreter","latex")
    ylabel("Stress $\sigma$ [MPa]", "Interpreter","latex")
    grid on
    xlim([-0.01,.4])
    ylim([0 350])
    %Extensometer Removed Location
    plot(S_Eng_Strain(IDX_S_E_Removal,1), ...
    S_Eng_Stress(IDX_S_E_Removal,1), "Marker","d","Color","b", ...
    "MarkerSize",6, "LineStyle","none")
    %Fracture Location
    plot(S_Eng_Strain(IDX_S_Fracture,1), ...
    S_Eng_Stress(IDX_S_Fracture,1), "Marker","x","Color","k", ...
    "MarkerSize",10, "LineStyle","none")
    plot(S_True_Strain(end,1), ...
    S_True_Stress(end,1), "Marker","x","Color","r", ...
    "MarkerSize",10, "LineStyle","none")
    % 0.2% Offset
    plot(S_Offset_Strain+0.002, S_Offset_Stress,"Color", "#7E2F8E", ...
    "LineStyle","-.", "LineWidth",1)
    %POI
    annotation('textarrow',[0.4 0.2],[0.35 0.45],'String','Yield',"Interpreter","latex")
    %LEGEND
    legend("Engineering Stress-Strain", "True Stress-Strain",...
     "Extensometer Removed","Engineering Fracture", "True Fracture", ...
     "0.2% Yield Offset", "Location", "Southeast")
%Plotting Aluminum
    %Stress and Strain
    figure(2)
    hold on
    plot(A_Eng_Strain, A_Eng_Stress, "Color","k","LineStyle","-")
    plot(A_True_Strain, A_True_Stress, "Color","r","LineStyle","--")
    title("6061 Aluminum Stress-Strain Curves")
    xlabel("Strain $\epsilon$  [mm/mm]","Interpreter","latex")
    ylabel("Stress $\sigma$ [MPa]", "Interpreter","latex")
    grid on
    xlim([-0.01,.1])
    ylim([0 350])
    %Extensometer Removed Location
    plot(A_Eng_Strain(IDX_A_E_Removal,1), ...
    A_Eng_Stress(IDX_A_E_Removal,1), "Marker","d","Color","b", ...
    "MarkerSize",6, "LineStyle","none")
    %Fracture Location
    plot(A_Eng_Strain(IDX_A_Fracture,1), ...
    A_Eng_Stress(IDX_A_Fracture,1), "Marker","x","Color","k", ...
    "MarkerSize",10, "LineStyle","none")
    plot(A_True_Strain(end,1), ...
    A_True_Stress(end,1), "Marker","x","Color","r", ...
    "MarkerSize",10, "LineStyle","none")
    % 0.2% Offset
    plot(A_Offset_Strain + 0.002, A_Offset_Stress,"Color", "#7E2F8E", ...
    "LineStyle","-.", "LineWidth",1)
    %POI
    annotation('textarrow',[0.4 0.28],[0.4 0.6],'String','Yield',"Interpreter","latex")
    %LEGEND
    legend("Engineering Stress-Strain", "True Stress-Strain",...
     "Extensometer Removed","Engineering Fracture", "True Fracture", ...
     "0.2% Yield Offset", "Location", "Southeast")
%Plotting DIC
    %Stress and Strain
    figure(3)
    hold on
    plot(DIC_Eng_Strain, DIC_Eng_Stress, "Color","k","LineStyle","-")
    plot(DIC_True_Strain, DIC_True_Stress, "Color","r","LineStyle","--")
    title("Digital Image Correlation Stress-Strain Curves")
    xlabel("Strain $\epsilon$  [mm/mm]","Interpreter","latex")
    ylabel("Stress $\sigma$ [MPa]", "Interpreter","latex")
    grid on
    xlim([-0.01,.1])
    ylim([0 600])
    %Extensometer Removed Location
    plot(DIC_Eng_Strain(IDX_DIC_E_Removal,1), ...
    DIC_Eng_Stress(IDX_DIC_E_Removal,1), "Marker","d","Color","b", ...
    "MarkerSize",6, "LineStyle","none")
    %Fracture Location
    plot(DIC_Eng_Strain(IDX_DIC_Fracture,1), ...
    DIC_Eng_Stress(IDX_DIC_Fracture,1), "Marker","x","Color","k", ...
    "MarkerSize",10, "LineStyle","none")
    plot(DIC_True_Strain(end,1), ...
    DIC_True_Stress(end,1), "Marker","x","Color","r", ...
    "MarkerSize",10, "LineStyle","none")
    % 0.2% Offset
    plot(DIC_Offset_Strain + 0.002, DIC_Offset_Stress,"Color", "#7E2F8E", ...
    "LineStyle","-.", "LineWidth",1)
    %POI
    annotation('textarrow',[0.4 0.25],[0.3 0.4],'String','Yield',"Interpreter","latex")
    %LEGEND
    legend("Engineering Stress-Strain", "True Stress-Strain",...
     "Extensometer Removed","Engineering Fracture", "True Fracture", ...
     "0.2% Yield Offset", "Location", "Southeast")
figure(4)
hold on
plot(ARC_LENGTH, VSG_STRAIN_YY, "Color","k","LineStyle","-")
title("Digital Image Correlation - YY VSG Strain")
xlabel("Arc length (pixels)","Interpreter","latex")
ylabel("Strain $\epsilon_{YY}$", "Interpreter","latex")
grid on

figure(4)
hold on
plot(ARC_LENGTH, VSG_STRAIN_XX, "Color","k","LineStyle","-")
title("Digital Image Correlation - YY VSG Strain")
xlabel("Arc length (pixels)","Interpreter","latex")
ylabel("Strain $\epsilon_{XX}$", "Interpreter","latex")
grid on

%Calibration 
figure(5)
plot(Calibration_Pos_Extensometer,Calibration_V)

%Plot Exporting
exportgraphics(figure(1), 'A_StressStrain.pdf', 'ContentType', 'vector','Resolution',300);
exportgraphics(figure(2), 'S_StressStrain.pdf', 'ContentType', 'vector','Resolution',300);
exportgraphics(figure(3), 'ADIC_StressStrain.pdf', 'ContentType', 'vector','Resolution',300);
exportgraphics(figure(4), 'ADIC_VSGStrain.pdf', 'ContentType', 'vector','Resolution',300);

%% Error Propagation
%Uncertainty values, all +-
    S_Load = 22.2+0.04*S_Force; %N
    A_Load = 22.2+0.04*A_Force; %N
    DIC_Load = 22.2+0.04*DIC_Force; %N
    CreepLoad = 2.22; %N
    CrossheadPos = 0.254; %mm
    S_Ext = 0.04*S_Extensometer; %mm
    A_Ext = 0.04*A_Extensometer; %mm
    DIC_Ext = 0.04*DIC_Extensometer; %mm
    Caliper = 0.127; %mm
%Steel
    S_Y_Err= (2.265-0.8723)*10^5/4;
    WS_stress = sqrt((S_Load./A0_S).^2+(-S_Force*Caliper./(A0_S*w_S_init)).^2+ ...
        (-S_Force*Caliper./(A0_S*t_S_init)).^2);
    WS_yield = WS_stress(21);
    WS_ult = max(WS_stress);
    WS_area = Caliper*sqrt((-t_S_fin/A0_S)^2+(-w_S_fin/A0_S)^2+(Afin_S/ ...
        (A0_S*t_S_init))^2+(Afin_S/(A0_S*w_S_init))^2);
    WS_TS = sqrt(CrossheadPos^2*(1/L_S(end))^2+Caliper^2*(1/GaugeL_S)^2);
%Aluminum
    A_Y_Err= (5.58-5.188)*10^4/4;
    WA_stress = sqrt((A_Load./A0_A).^2+(-A_Force*Caliper./(A0_A*w_A_init)).^2+ ...
        (-A_Force*Caliper./(A0_A*t_A_init)).^2);
    WA_yield = WA_stress(66);
    WA_ult = max(WA_stress);
    WA_area = Caliper*sqrt((-t_A_fin/A0_A)^2+(-w_A_fin/A0_A)^2+(Afin_A/ ...
        (A0_A*t_A_init))^2+(Afin_A/(A0_A*w_A_init))^2);
    WA_TS = sqrt(CrossheadPos^2*(1/L_A(end))^2+Caliper^2*(1/GaugeL_A)^2);
%DIC
    DIC_Y_Err = (3.577-3.075)*10^5/4;
    WDIC_stress = sqrt((DIC_Load./A0_DIC).^2+(-DIC_Force*Caliper./ ...
        (A0_DIC*w_DIC_init)).^2+(-DIC_Force*Caliper./(A0_DIC*t_DIC_init)).^2);
    WDIC_yield = WDIC_stress(38);
    WDIC_ult = max(WDIC_stress);
    WDIC_area = Caliper*sqrt((-t_DIC_fin/A0_DIC)^2+(-w_DIC_fin/A0_DIC)^2+(Afin_DIC/ ...
        (A0_DIC*t_DIC_init))^2+(Afin_DIC/(A0_DIC*w_DIC_init))^2);
    WDIC_TS = sqrt(CrossheadPos^2*(1/L_DIC(end))^2+Caliper^2*(1/GaugeL_DIC)^2);