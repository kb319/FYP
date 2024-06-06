close all;
clear all;
clc;

COP_data = readtable('HP_COP.csv');
COP_data = table2array(COP_data(:,:));

capacity_data = readtable('HP_capacity.csv');
capacity_data = table2array(capacity_data(:,:));

row_ambient_temp = [-20 -15 -10 -7 2 7 12 15 20];
column_flow_temp = [25 35 40 45 50 55 60];


[grid_flow, grid_ambient] = meshgrid( column_flow_temp, row_ambient_temp);

p = polyfitn([grid_flow(:), grid_ambient(:)], COP_data(:), 2); % Fit a quadratic polynomial NEEDS TO BE AT LEAST CUBIC LATER
Z_fit = polyvaln(p, [grid_flow(:), grid_ambient(:)]);
Z_fit = reshape(Z_fit, size(grid_ambient));


% First figure
figure;
surf(grid_ambient, grid_flow, COP_data);
ax = gca;
ax.LineWidth = 2;
ax.FontSize = 14; 

xlabel('Ambient Temperature (°C)', 'FontSize', 14, 'FontWeight', 'bold', 'Position', [mean(get(gca, 'xlim'))-3, min(get(gca, 'ylim')) - 0.5, min(get(gca, 'zlim'))-1.4], 'Rotation', 16);
ylabel('Flow Temperature (°C)', 'FontSize', 14, 'FontWeight', 'bold', 'Position', [min(get(gca, 'xlim')) - 0.5, mean(get(gca, 'ylim'))+1, min(get(gca, 'zlim'))-1.4], 'Rotation', -25);
zlabel('Coefficient of Performance', 'FontSize', 14, 'FontWeight', 'bold', 'Rotation', 90);

% Second figure
figure;
surf(grid_ambient, grid_flow, capacity_data);
ax = gca;
ax.LineWidth = 2; 
ax.FontSize = 14; 

xlabel('Ambient Temperature (°C)', 'FontSize', 14, 'FontWeight', 'bold', 'Position', [mean(get(gca, 'xlim'))-3, min(get(gca, 'ylim')) - 0.5, min(get(gca, 'zlim'))-0.4], 'Rotation', 16);
ylabel('Flow Temperature (°C)', 'FontSize', 14, 'FontWeight', 'bold', 'Position', [min(get(gca, 'xlim')) - 0.5, mean(get(gca, 'ylim'))+1, min(get(gca, 'zlim'))-0.4], 'Rotation', -25);
zlabel('Capacity (kW)', 'FontSize', 14, 'FontWeight', 'bold', 'Rotation', 90);

