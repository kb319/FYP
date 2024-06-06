clear all;
close all;
clc; 
tic;
 
day_of_year = 1; % KEY 

day_of_year = (day_of_year -1)*2*24;

sim_days = 364;
Ts = 0.5; % 1/amount of samples per hour
n = 24/Ts*sim_days; % sim length in points
horz_length_days = 2;
%N = horz_length_days/sim_days*n; % horizon, so 48 is 24hrs

C_a = 0.23; % adjust so loss is around 1 degree an hour with no input at day peak
U_a = 0.16;
C_h = 0.20;
U_h = 0.14;


pricedata = readtable('octop.csv');
p_buy_price = table2array(pricedata(:,1))';


tspan = linspace(0,sim_days*24*60*60,n); %1 minute intervals

outdoor_temp = readtable('weather_thirtymin.csv');
outdoor_temp = table2array(outdoor_temp(1,:));
%% mapping for flow and indoor temp to get COP

COP_data = readtable('HP_COP.csv');
COP_data = table2array(COP_data(:,:));

capacity_data = readtable('HP_capacity.csv');
capacity_data = table2array(capacity_data(:,:));

row_ambient_temp = [-20 -15 -10 -7 2 7 12 15 20];
column_flow_temp = [25 35 40 45 50 55 60];

[grid_flow, grid_ambient] = meshgrid( column_flow_temp, row_ambient_temp);

p = polyfitn([grid_flow(:), grid_ambient(:)], COP_data(:), 2); 

%% home RC model and discretising
T_s = Ts;

A_home = [-(U_a+U_h)/C_a U_h/C_a; U_h/C_h -U_h/C_h];
B_home = [U_a/C_a 0; 0 1/C_h];

syshome = idss(A_home,B_home,[],[],[],[],0);
syshome_disc = c2d(syshome,T_s,'zoh');

initial_indoor_temp = 21;
initial_flow_temp =44;

A_home_disc = syshome_disc.A;
B_home_disc = syshome_disc.B;
A = A_home_disc;
B = B_home_disc;
%initial state
x1_o = [initial_indoor_temp initial_flow_temp];
    T_in_plot = [x1_o(1)];
    T_flow_plot = [x1_o(2)];
    Q_HP_prev = 6;
    Q_HP_plot = [];
    COP_plot = [];
    psi_plot = [];
    T_o_plot = [];
for i = 1+day_of_year:n+day_of_year
    if x1_o(1) >= 22
        Q_HP = 0;
    elseif x1_o(1) <= 20
        Q_HP = 6;

    else 
        Q_HP = Q_HP_prev;
    end

    if Q_HP_prev == 0 && Q_HP == 6
        COP_new = polyvaln(p,[x1_o(2),outdoor_temp(i)])*0.95;
    else 
        COP_new = polyvaln(p,[x1_o(2),outdoor_temp(i)]);
    end

    psi_new = Q_HP/COP_new;

    COP_plot = [COP_plot COP_new];
    psi_plot = [psi_plot psi_new];


    y = A*x1_o'+B*[outdoor_temp(i) Q_HP]';
    x1_o = y';


    Q_HP_plot = [Q_HP_plot Q_HP];
    T_in_plot = [T_in_plot x1_o(1)];
    T_flow_plot = [T_flow_plot x1_o(2)];
    T_o_new = outdoor_temp(i);
   T_o_plot = [T_o_plot T_o_new];
    
    Q_HP_prev = Q_HP; 
end

toc;

fprintf('total psi = %f\n',(sum(psi_plot)));

linesize = 0.1;

figure;
set(gca, 'FontSize', 12, 'LineWidth', linesize); 
hold on;
grid on;
%correct_x_gaps = linspace(0,sim_days*24,n);
plot(1:size(T_o_plot,2),T_o_plot,'LineWidth',linesize);
plot(1:size(T_in_plot,2),T_in_plot,'LineWidth',linesize);
plot(1:size(T_flow_plot,2),T_flow_plot,'LineWidth',linesize);
plot(1:size(COP_plot,2),COP_plot,'LineWidth',linesize);
yline(22,'--');
yline(20,'--');
% making the temp and input datasets interpolatable
xlabel("Data Index (half-hourly intervals)");
ylabel("Temperature (°C)")
yyaxis right;
ylabel("Power (kW)");
stairs(1:size(Q_HP_plot,2),Q_HP_plot,'LineWidth',linesize);
yline(0,'--','color','g');
yline(6,'--','color','g');
ylim([0 50])
legend("Outdoor Temperature","Indoor Temperature","Flow Temperature","Coefficient of Performance","","","Heat Pump Thermal Output","","",'FontSize', 12, 'Location', 'north')

p_buy_price_calc = p_buy_price(1+day_of_year:size(psi_plot,2)+day_of_year);

tot_buy_calc = sum(p_buy_price_calc.*psi_plot*Ts);
tot_pounds = (tot_buy_calc)/100;
fprintf('Total Compressor Power = %f\n',(sum(psi_plot)));
fprintf('Total Cost = £%f\n',tot_pounds);