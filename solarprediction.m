clear all;
close all;
clc;

%% extracting data, normalising to scale of 1
pv_real = readtable('P_pv_thirtymin.csv');
pv_real = table2array(pv_real(1,:));
pv_real_norm = pv_real./max(pv_real(:));

pv_atmos = readtable('P_pv_atmos_thirtymin.csv');
pv_atmos = table2array(pv_atmos(1,:));
pv_atmos_norm = pv_atmos./max(pv_atmos(:));

cloud = readtable('cloud_thirtymin.csv');
cloud = table2array(cloud(1,:));

%% making non-light hours cloud data zero
pv_atmos_saturated = pv_atmos~=0;
cloud_adj = cloud.*pv_atmos_saturated;

%% fitting a quadratic to data relating okta 1-cloud cover to solar irradiance

x_data = [8 7 6 5 4 3 2 1 0];
x_data = x_data./max(x_data(:));
y_data = [200 300 410 510 610 680 730 750 780];
y_data = y_data./max(y_data(:));
%89y_data = flip(y_data,2);
% Fit a quadratic curve (2nd degree polynomial) to the data
p_coefficients = polyfit(x_data, y_data, 2);
% Generate x values for plotting the fitted curve
x_fit = linspace(min(x_data), max(x_data), 100); 

y_fit = polyval(p_coefficients, x_fit);
%y_fit = flip(y_fit,2);
okta_vec = 0:8;
okta_vec_fit = linspace(0,8,100);
% Plot original datafigure;
figure;
set(gca, 'FontSize', 12, 'LineWidth', 1.2); % Set axes properties
hold on;
grid on;
plot(okta_vec, y_data, 'bo', 'DisplayName', 'Original Data','LineWidth',1.2);
plot(okta_vec_fit, y_fit, 'r-', 'DisplayName', 'Quadratic Fit','LineWidth',1.2);
xlabel("Cloud Cover (okta)");
ylabel("Normalised Solar Irradiance (relative)")
legend("Measured Datapoints","Quadratic Fit")

%% testing out solar prediction

pv_pred = pv_atmos_norm.*polyval(p_coefficients,cloud_adj).*max(pv_real);
figure;
set(gca, 'FontSize', 14, 'LineWidth', 2); % Set axes properties

% Define the length
length = 1:17000;

subplot(2, 1, 1);
set(gca, 'FontSize', 12, 'LineWidth', 1.2);
hold on;
grid on;
ylim([0 5]);
plot(length, pv_real(length), 'LineWidth', 1.2);
plot(length, pv_atmos_norm(length), 'LineWidth', 1.2);
%plot(length, pv_pred(length), 'LineWidth', 1.2);
ylabel("Solar Array Output (kW)");
yline(3.65);
xlabel("Data Index (Half-Hourly Interval)");
legend("Real Solar Array Output", "Atmospheric Solar Irradiation", "Solar Array Output Prediction","location",'north');
hold off;

subplot(2, 1, 2);
set(gca, 'FontSize', 12, 'LineWidth', 1.2);
hold on;
grid on;
plot(length, cloud(length), 'LineWidth', 1.2);
ylim([0 1]);
ylabel("Scaled Cloud Cover Level");
xlabel("Data Index (Half-Hourly Interval)");
legend("True Cloud Cover","location",'northwest');
hold off;

% Calculate mean squared errors
err_exact = immse(pv_pred, pv_real);
% Display errors
disp(['Mean Squared Error: ', num2str(err_exact)]);



