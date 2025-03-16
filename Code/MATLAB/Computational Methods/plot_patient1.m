%PLOT_PATIENT1 Saves the daily average inflation
patient_data = readmatrix('matlab-basics-data\inflammation-01.csv');

%Plot the daily average for all patients
figure
plot(mean(patient_data, 1));
title('Daily Average Inflammation')
xlabel('Day of Trial')
ylabel('Inflammation')

% Save figure
print('matlab-basics-data/daily_avg.png', '-dpng');
