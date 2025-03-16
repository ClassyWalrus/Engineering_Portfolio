close all;

replace = 0;

% create a set of alphas and betas over range of values
linSpaceAlphas = 0:.2:25;
linSpaceBetas = 0:.5:50;
% create grid
[alphas,betas] = meshgrid(linSpaceAlphas,linSpaceBetas);
% make the grids into arrays (instead of matrices)
alphas = alphas(:); betas = betas(:);
% get rid of "infeasible" combinations (where point C is further than point
% point B)
infeasible =  (alphas > betas);
alphas(infeasible) =[];
betas(infeasible) = [];

% calculate the delta V for the Hohmann and the bi-elliptic Hohmann;
% Us the "alphaBoundHohmMoreEfficient", "alphaBoundBiMoreEfficient",
% "betas", and "alphas" variables in order to do so
       % HINT: it will be helpful to use the "element-wise" operation (e.g. ".*",
       % "./", ".^" within these lines to figure out based on the alphas and betas
delta_v_Hohmann = 1./(sqrt(alphas))-((sqrt(2).*(1-alphas)))./(sqrt(alphas.*(1+alphas)))-1;
delta_v_biHohmann = sqrt((2.*(alphas+betas))./(alphas.*betas))-((1+sqrt(alphas))./sqrt(alphas))-sqrt(2./(betas.*(1+betas))).*(1-betas);
% 1/(sqrt(alphas))-((sqrt(2)*(1-alphas)))/(sqrt(alphas.*(1+alphas)));
% determine the advantage of performing the Hohmann transfer over the
% bielliptic Hohmann (or in some cases, disadvantage)
% by calculating the additional delta-V needed to perform the bielliptic
% Hohmann transfer (negative values where bielliptic Hohmann is more
% efficient)
delta_v_HohmannAdvantage = delta_v_biHohmann-delta_v_Hohmann;


% set up colors for scatter plot
lightBlue =[0 255 255]/255;
lightYellow = [255 255 0]/255;
colors_v_Hohmann = repmat(lightBlue, length(alphas),1);
colors_v_biHohmann = repmat(lightYellow, length(alphas),1);



% Two scatter plots to show...
figure;
% dont worry about the fact that I plotted -delta_v values (because I set
% the tick labels latter)
scatter3(alphas, betas, -delta_v_Hohmann, [],colors_v_Hohmann,'filled');
hold on;
scatter3(alphas, betas, -delta_v_biHohmann, [],colors_v_biHohmann,'filled');
zt = get(gca, 'ZTick');
set(gca, 'ZTick',zt, 'ZTickLabel',-zt)
xlabel('\alpha'); ylabel ('\beta'); zlabel ('\Delta V')

% Set up value for when the Hohmann transfer is more efficient (regardless of B point of bielliptic)
alphaBoundHohmMoreEfficient = 11.94;
% Set up value for when bielliptic transfer is more efficient (regardless of B point of bielliptic)
alphaBoundBiMoreEfficient = 15.58;

% range of beta values associated with alpha values above
betasHohm = linSpaceBetas(linSpaceBetas>alphaBoundHohmMoreEfficient);
betasBi = linSpaceBetas(linSpaceBetas>alphaBoundBiMoreEfficient);

% Similar to other calculation of deltaVs above, calculate the delta V for the Hohmann and the bi-elliptic Hohmann;
% Use the "alphaBoundHohmMoreEfficient", "alphaBoundBiMoreEfficient",
% "betasHohm", and "betasBi" variable in order to do so
       % HINT: it will be helpful to use the "element-wise" operation (e.g. ".*",
       % "./", ".^" within these lines to figure out based on the alphas and betas
delta_v_Hohmann_Bound = 1./(sqrt(alphaBoundHohmMoreEfficient))-((sqrt(2).*(1-alphaBoundHohmMoreEfficient)))./(sqrt(alphaBoundHohmMoreEfficient.*(1+alphaBoundHohmMoreEfficient)))-1;
delta_v_biHohmann_Bound = sqrt((2.*(alphaBoundBiMoreEfficient+betasBi))./(alphaBoundBiMoreEfficient.*betasBi))-((1+sqrt(alphaBoundBiMoreEfficient))./sqrt(alphaBoundBiMoreEfficient))-sqrt(2./(betasBi.*(1+betasBi))).*(1-betasBi);

% black lines plotted at alpha boundaries
scatter3(alphaBoundHohmMoreEfficient*ones(length(betasHohm),1), betasHohm, -delta_v_Hohmann_Bound, [],'k','filled');
scatter3(alphaBoundBiMoreEfficient*ones(length(betasBi),1), betasBi, -delta_v_biHohmann_Bound, [],'k','filled');

% text labels
textLabel1 = 'Hohmann More Efficient';
textLabel2 = 'Bielliptical More Efficient';
text(alphaBoundHohmMoreEfficient, mean(betasHohm), textLabel1,'HorizontalAlignment','right');
text(alphaBoundBiMoreEfficient, mean(betasBi), textLabel2,'HorizontalAlignment','left');
view(0,90);

figure;
scatter3(alphas, betas, delta_v_HohmannAdvantage,[],log(abs(delta_v_HohmannAdvantage).^-100),'filled');
hold on;
scatter3(alphas(abs(delta_v_HohmannAdvantage)<1e-4), betas(abs(delta_v_HohmannAdvantage)<1e-4), delta_v_HohmannAdvantage(abs(delta_v_HohmannAdvantage)<1e-4),[],'k','filled');
xlabel('\alpha'); ylabel ('\beta'); zlabel ('\Delta V difference')


time_HohmannAdvantage = ((1+alphas).^(3/2))./((1+betas).^(3/2)+(betas+alphas).^(3/2));
figure;
colormap(summer);
scatter3(alphas(delta_v_HohmannAdvantage<0), betas(delta_v_HohmannAdvantage<0), time_HohmannAdvantage(delta_v_HohmannAdvantage<0),[],time_HohmannAdvantage(delta_v_HohmannAdvantage<0),'filled');
hold on;
scatter3(alphas(abs(delta_v_HohmannAdvantage)<1e-4), betas(abs(delta_v_HohmannAdvantage)<1e-4), time_HohmannAdvantage(abs(delta_v_HohmannAdvantage)<1e-4),[],'k','filled');
xlabel('\alpha'); ylabel ('\beta'); zlabel ('time(bi)/time(H)')
colorbar;
cb = colorbar;
cb.Label.String = 'tbi/th';