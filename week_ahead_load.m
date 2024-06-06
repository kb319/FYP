clear all;
clc; 

load = readtable('load_thirtymin.csv');

load = table2array(load(1,:));

week = 48*7;
load_past = load(1:end-week);
load_shift = load(1+week:end);

figure;
hold on;

stairs(1:17184,load_past);
stairs(1:17184,load_shift);

MSE = immse(load_past,load_shift);

figure;
set(gca, 'FontSize', 12, 'LineWidth', 1.2); % Set axes properties
hold on;
grid on;
stairs(1:17184, load_past,'LineWidth',1.2);
stairs(1:17184, load_shift,'--','LineWidth',1.2);
xlabel("Data Index (half-hourly intervals)");
ylabel("Power (kW)")
legend("Previous Week Shifted Electrical Load","Electrical Load","Location","best")

