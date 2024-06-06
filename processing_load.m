clear all;
clc; 

load = readtable('home68_filled_one_year_hourly.csv');

load = table2array(load(:,1))';


load_app = [load load(end)];

load_avg = [];
for i = 1:size(load,2)
    val = (load_app(i)+load_app(i+1))/2;
    load_avg = [load_avg val];
end

load_30 = vertcat(load,load_avg);

load_30 = reshape(load_30,[1,17520]);
load_30 = load_30./1000;
figure;
plot(1:8760,load);

figure;
plot(1:17520,load_30);

writematrix(load_30,'load_thirtymin.csv');