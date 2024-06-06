clear all;
clc; 

P_pv = readtable('pv_imperial.csv');

P_pv = table2array(P_pv(:,1))';

plot(1:10,P_pv(8:17));

P_pv_app = [P_pv P_pv(end)];

P_avg = [];
for i = 1:size(P_pv,2)
    val = (P_pv_app(i)+P_pv_app(i+1))/2;
    P_avg = [P_avg val];
end

P_pv_30 = vertcat(P_pv,P_avg);

P_pv_30 = reshape(P_pv_30,[1,17520]);

figure;
plot(1:8760,P_pv);

figure;
plot(1:17520,P_pv_30);

writematrix(P_pv_30,'P_pv_thirtymin.csv');