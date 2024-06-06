clear all;
close all;
clc;
tic;

%% Battery characteristics in kW
P_buy_max = 5;
P_sell_max = 5;
P_charge_max = 3.68;
P_discharge_max = 5;
Cmin = 0.1;
Cmax = 0.9;
eta_ch = 0.95;
eta_dis = 0.95;
bat_capacity = 10;  

%% Setting timeframes

day_of_year = 180; % KEY 

day_of_year = (day_of_year -1)*2*24;

sim_days = 2;

Ts = 0.5; % 1/amount of samples per hour
n = 24/Ts*sim_days; % sim length in points
horz_length_days = 2;
global N;
N = horz_length_days*24/Ts; % horizon, so 48 is 24hrs

%% Creating pv and load data, loading buy and sell

P_load = readtable('load_thirtymin.csv');
P_load = table2array(P_load(1,:));


t_samp = linspace(0,24*sim_days-0.5,n);

P_pv = readtable('P_pv_thirtymin.csv');
P_pv = table2array(P_pv(1,:));
P_pv = P_pv(1:size(P_load,2));

pricedata = readtable('octop.csv');
global p_buy_price;
global p_sell_price;
p_buy_price = table2array(pricedata(:,1))';
p_sell_price =  table2array(pricedata(:,2))';

figure;
hold on;

set(gca, 'FontSize', 12, 'LineWidth', 1.2); % Set font size and line width for the current axes
plot(1:size(P_pv, 2), P_load, 'LineWidth', 1.5); % Plot with thicker lines
plot(1:size(P_pv, 2), P_pv, 'LineWidth', 1.5); 

legend('P_{load}', 'P_{pv}', 'FontSize', 12, 'Location', 'best'); % Add a legend with font size
xlabel('Time (hours)', 'FontSize', 12, 'FontWeight', 'bold'); % Add x label with font size and weight
ylabel('Normalised Power', 'FontSize', 12, 'FontWeight', 'bold'); % Add y label with font size and weight

% Apply general settings for the entire figure
set(findall(gcf, '-property', 'FontSize'), 'FontSize', 12); % Ensure all text has the specified font size
set(findall(gcf, '-property', 'FontWeight'), 'FontWeight', 'bold'); % Ensure all text is bold

%% Battery dynamics + initial and target conditions
A = 0.99;
B = [0 0 Ts*eta_ch -Ts*eta_dis];

initial_charge = 0.4;
y = initial_charge*bat_capacity;
target_charge = 0.9;
yend = target_charge*bat_capacity;

%% Controller
control_mode = 1; % 1 Predictive Control, 2 self-powered, 3 PC with PV prediction 
% 4 PC with PV prediction and constant load prediction
switch control_mode
    case 1
        global x_size;
        x_size = size(A,1);
        u_size = size(B,2);
        
        % Doing Aeq 
        
        Aeq_AI = horzcat(kron(eye(N),A),zeros(N,x_size)) - horzcat(zeros(N,x_size),kron(eye(N),eye(x_size)));
        Aeq_B = kron(eye(N),B);
        Aeq_balance = horzcat(zeros(N,N+1),kron(eye(N),[1 -1 -1 1]));
        Aeq_base = [Aeq_AI Aeq_B;Aeq_balance];
        
        % input upper and lower bounds
        
        lb_u = [0 0 0 0];
        ub_u = [P_buy_max P_sell_max P_charge_max P_discharge_max];
        lb_x = Cmin*bat_capacity;
        ub_x = Cmax*bat_capacity;
        lb = [kron(ones(1,N+1),lb_x) kron(ones(1,N),lb_u)];
        ub = [kron(ones(1,N+1),ub_x) kron(ones(1,N),ub_u)];
        
        charge_state = [];
        buy_pow = [];
        sell_pow = [];
        ch_pow = [];
        dis_pow = [];
        
        %running mpc  

        initial_guess = ones(1,(N+1)*x_size+N*u_size);
        global size_decision;
        size_decision = size(initial_guess,2);
        global i;
        options = optimoptions("fmincon","Algorithm","interior-point","EnableFeasibilityMode",false,"MaxFunctionEvaluations",3000000,"ConstraintTolerance",1e-6,"MaxIterations",1e7);

        for i = 1+day_of_year:n+day_of_year
        
            Beq_base = [zeros(N*x_size,1); transpose(P_load(i:i+N-1) - P_pv(i:i+N-1))];
            Aeq = vertcat([eye(x_size),zeros(1,N*x_size+N*u_size)],[zeros(1,N),1,zeros(1,N*4)],Aeq_base);
            Beq = vertcat(y, yend, Beq_base);
        
            % doing  cost function
            
            sols = fmincon(@(d) cost_function(d),initial_guess,[],[],Aeq,Beq,lb,ub,[],options);
            initial_guess = sols;
         
            u = sols(N+2:N+5);
            x = sols(1);
            y = A*x+B*u';
        
            % analysis
        
            charge_state = [charge_state x];
            buy_pow = [buy_pow u(1)];
            sell_pow = [sell_pow u(2)];
            ch_pow = [ch_pow u(3)];
            dis_pow = [dis_pow u(4)];
        end

end

%% Calculating and displaying total cost 

p_buy_price_calc = p_buy_price(1+day_of_year:size(buy_pow,2)+day_of_year);
p_sell_price_calc = p_sell_price(1+day_of_year:size(buy_pow,2)+day_of_year);

tot_buy_calc = sum(p_buy_price_calc.*buy_pow*Ts);
tot_sell_calc = sum(p_sell_price_calc.*sell_pow*Ts);
tot_pounds = (tot_buy_calc-tot_sell_calc)/100;
fprintf('Total cost = £%f\n',tot_pounds);

%% Plotting results
figure;
set(gca, 'FontSize', 14, 'LineWidth', 1.2); % Set axes properties
plotlength = linspace(0, n * Ts, n);

subplot(3,1,1)
set(gca, 'FontSize', 14, 'LineWidth', 1.2); % Set axes properties
plotlength = linspace(0, n * Ts, n);
hold on;
plot(plotlength, buy_pow, 'LineWidth', 1.5); % Thicker plot lines
plot(plotlength, sell_pow, 'LineWidth', 1.5); 
plot(plotlength, ch_pow - dis_pow, 'LineWidth', 1.5); 
xlabel('Simulation Hour', 'FontSize', 12);
ylabel('Power (kW)', 'FontSize', 12);
legend('Buy Power', 'Sell Power', 'Charge/Discharge Power', 'FontSize', 12, 'Location', 'best');
set(gca, 'FontSize', 15, 'LineWidth', 1.2);

subplot(3,1,2)
set(gca, 'FontSize', 14, 'LineWidth', 1.2); % Set axes properties
plotlength = linspace(0, n * Ts, n);
hold on;
plot(plotlength, charge_state, 'LineWidth', 1.2);
yline(1, '--', 'LineWidth', 1.2);
yline(9, '--', 'LineWidth', 1.2);
xlabel('Simulation Hour', 'FontSize', 12);
ylabel('Charge Level (kWh)', 'FontSize', 12);
legend('Battery Charge', 'Minimum Charge', 'Maximum Charge', 'FontSize', 12, 'Location', 'best');
set(gca, 'FontSize', 15, 'LineWidth', 1.2);

subplot(3,1,3)
set(gca, 'FontSize', 14, 'LineWidth', 1.2); % Set axes properties
plotlength = linspace(0, n * Ts, n);
hold on;
plot(plotlength, p_buy_price(1+day_of_year:n+day_of_year), 'LineWidth', 1.2);
plot(plotlength, p_sell_price(1+day_of_year:n+day_of_year), 'LineWidth', 1.2);
xlabel('Simulation Hour', 'FontSize', 12);
ylabel('Prices (p/kWh)', 'FontSize', 12);
legend('Buy Price', 'Sell Price', 'FontSize', 12, 'Location', 'best');
set(gca, 'FontSize', 15, 'LineWidth', 1.2);


toc;

function objective = cost_function(d)
   global p_buy_price;
   global p_sell_price;
   global size_decision;
   global N;
   global i;
   global x_size;
   f = [];
        for j = 0:N-1
            new = [p_buy_price(i+j) -p_sell_price(i+j) 0 0];
            f = [f new]; 
        end
   f = horzcat(zeros(1,x_size*(N+1)),f); % putting zeros to skip the state part of the decision matrix variables
           

   objective = sum(f.*d(1:1:size_decision));
end
