clear all;
clc; 

P_cloud = readtable('cloud_imperial.csv');

P_cloud = table2array(P_cloud(:,1))';

plot(1:10,P_cloud(8:17));

P_cloud_app = [P_cloud P_cloud(end)];

P_avg = [];
for i = 1:size(P_cloud,2)
    val = (P_cloud_app(i)+P_cloud_app(i+1))/2;
    P_avg = [P_avg val];
end

P_cloud_30 = vertcat(P_cloud,P_avg);

P_cloud_30 = reshape(P_cloud_30,[1,17520]);

figure;
plot(1:8760,P_cloud);

figure;
plot(1:17520,P_cloud_30);

writematrix(P_cloud_30,'cloud_thirtymin.csv');