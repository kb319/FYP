clear all;
clc; 

P_temp = readtable('outdoortemp_imperial.csv');

P_temp = table2array(P_temp(:,1))';

plot(1:10,P_temp(8:17));

P_temp_app = [P_temp P_temp(end)];

P_avg = [];
for i = 1:size(P_temp,2)
    val = (P_temp_app(i)+P_temp_app(i+1))/2;
    P_avg = [P_avg val];
end

P_temp_30 = vertcat(P_temp,P_avg);

P_temp_30 = reshape(P_temp_30,[1,17520]);

figure;
plot(1:8760,P_temp);

figure;
plot(1:17520,P_temp_30);

writematrix(P_temp_30,'weather_thirtymin.csv');