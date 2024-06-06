close all;
clear all;
clc;
tic;
%% loading price data


pricedata = readtable('octop.csv');
global p_buy_price;
p_buy_price = table2array(pricedata(:,1))';

%% making COP mapping 
% 
global p;
COP_data = readtable('HP_COP.csv');
COP_data = table2array(COP_data(:,:));

capacity_data = readtable('HP_capacity.csv');
capacity_data = table2array(capacity_data(:,:));

row_ambient_temp = [-20 -15 -10 -7 2 7 12 15 20];
column_flow_temp = [25 35 40 45 50 55 60];


[grid_flow, grid_ambient] = meshgrid( column_flow_temp, row_ambient_temp);

p = polyfitn([grid_flow(:), grid_ambient(:)], COP_data(:), 2); % Fit a quadratic polynomial NEEDS TO BE AT LEAST CUBIC LATER
% Z_fit = polyvaln(p, [grid_flow(:), grid_ambient(:)]);
% Z_fit = reshape(Z_fit, size(grid_ambient));
% 
% 
% figure;
% surf(grid_ambient, grid_flow, COP_data);
% figure;
% surf(grid_ambient, grid_flow, capacity_data);



%% MPC heating

day_of_year = 10; % KEY 

day_of_year = (day_of_year -1)*2*24;


sim_days =2;
Ts = 0.5; % 1/amount of samples per hour
n = 24/Ts*sim_days; % sim length in points
horz_length_days = 1;
%N = horz_length_days/sim_days*n; % horizon, so 48 is 24hrs



tspan = linspace(0,sim_days*24*60*60,n); %1 minute intervals

outdoor_temp = readtable('weather_thirtymin.csv');
outdoor_temp = table2array(outdoor_temp(1,:));

C_a = 0.23; % adjust so loss is around 1 degree an hour with no input at day peak
U_a = 0.16;
C_h = 0.20;
U_h = 0.14;


%% home RC model and discretising
T_s = 0.5;

A_home = [-(U_a+U_h)/C_a U_h/C_a; U_h/C_h -U_h/C_h];
B_home = [U_a/C_a 0; 0 1/C_h];

syshome = idss(A_home,B_home,[],[],[],[],0);
syshome_disc = c2d(syshome,T_s,'zoh');

global target_indoor_temp;

initial_indoor_temp = 21;
target_indoor_temp = 21;

initial_flow_temp = 44.96-1.14*(outdoor_temp(1+day_of_year));
target_flow_temp = 50;



A_home_disc = syshome_disc.A;
B_home_disc = syshome_disc.B;
A = A_home_disc;
B = B_home_disc;
%initial state
x1_o = [initial_indoor_temp initial_flow_temp];

% doing MPC

global x_size
x_size = size(A,1);
u_size = size(B,2);

global N;
N = 20;
y = [initial_indoor_temp, initial_flow_temp]';
yend = [target_indoor_temp, target_flow_temp]'; 
Q_a = [];



% system dynamic equality constraints
Aeq_AI = horzcat(kron(eye(N),A),zeros(x_size*N,x_size)) - horzcat(zeros(N*x_size,x_size),kron(eye(N),eye(x_size)));
Aeq_B = kron(eye(N),B);
Aeq_base = [Aeq_AI Aeq_B];
Beq_base = [zeros(N*x_size,1)];

% setting outdoor temperature inputs as equality constraints
Aeq_temps = kron(eye(N),[1 zeros(1,u_size-1)]);
Aeq_temps = horzcat(zeros(N,x_size*(N+1)),Aeq_temps);
Aeq = vertcat([eye(x_size),zeros(x_size,N*x_size+(N)*u_size)],[zeros(x_size,(N)*x_size),eye(x_size),zeros(x_size,u_size*(N))],Aeq_base,Aeq_temps); % adding initial conditions
Aeq = horzcat(Aeq,zeros(size(Aeq,1),N));

%(Beq_temps in the main for loop because it shifts)

% Q inputs must be greater than zero and less than six

A_Q_zero = -kron(eye(N),[0 1 zeros(1,u_size-2)]);
A_Q_fit_zero = horzcat(zeros(N,x_size*(N+1)),A_Q_zero);
B_Q_zero = zeros(N,1);
 
A_Q_six = kron(eye(N),[0 1 zeros(1,u_size-2)]);
A_Q_fit_six = horzcat(zeros(N,x_size*(N+1)),A_Q_six);
B_Q_six = 6*ones(N,1);

A_Q = vertcat(A_Q_fit_zero,A_Q_fit_six);
B_Q = vertcat(B_Q_zero,B_Q_six);

% adding zeros for psi

A_Q = horzcat(A_Q,zeros(size(A_Q,1),N));
global size_decision
size_decision = size(Aeq,2);

% cost function

global initial_guess;
initial_guess = zeros(1,size_decision);
initial_guess(1:2*(N+1)) = kron(ones(1,N+1),[initial_indoor_temp,initial_flow_temp]);
initial_guess(2*(N+1)+1:size_decision-N) = kron(ones(1,N),[4,20]);
initial_guess(size_decision-N+1:end) = 2*ones(1,N);

options = optimoptions("fmincon","Algorithm","interior-point","EnableFeasibilityMode",false,"MaxFunctionEvaluations",300000,"ConstraintTolerance",1e-4,"MaxIterations",1e7);
T_in_plot = [initial_indoor_temp];
T_flow_plot = [initial_flow_temp];
Q_out = [];
Cop_plot = [];
psi_plot = [];
T_o_plot = [];
global i;
for i = 1+day_of_year:n+day_of_year

   i = i
   Beq_temps = outdoor_temp(i:i+N-1)';
  
   Beq = vertcat(y, yend, Beq_base, Beq_temps);
   % layout of decision variables is Indoor temps (X), then outdoor temp
   % and heat output alternating, then the compressor speeds
    sols = fmincon(@(d) cost_function(d),initial_guess,[],[],Aeq,Beq,[],[],@psi_to_q,options);
   %sols = fmincon(@(d) cost_function(d),initial_guess,A_Q,B_Q,Aeq,Beq,[],[],@psi_to_q,options);
   x1_o = sols(1:2)';             
   Q_HP_new = sols(2*(N+1)+2);
   T_o_sol = sols(2*(N+1)+1);
   
   initial_guess = sols;
   if Q_HP_new < 0
       Q_HP_new = 0;
       sols(size_decision-N+1) = 0;
   elseif  Q_HP_new > 6
       Q_HP_new = 6;
       sols(size_decision-N+1) = 6/polyvaln(p,[x1_o(2),outdoor_temp(i)]);
   end


   y = A*x1_o+B*[outdoor_temp(i) Q_HP_new]';

   troub_t_in = sols(1:2:2*(N+1)-1);
   troub_t_flow = sols(2:2:2*(N+1));
   troub_t_out = sols(2*(N+1)+1:2:size_decision-N-1);
   troub_q_out = sols(2*(N+1)+2:2:size_decision-N);
   Cop_plot = [Cop_plot polyvaln(p,[x1_o(2),outdoor_temp(i)])];
   
   psi_plot = [psi_plot sols(size_decision-N+1)];
   T_o_new = outdoor_temp(i);
   T_o_plot = [T_o_plot T_o_new];
   T_in_plot = [T_in_plot y(1)]; % y-target as the cost is based around it, so this gives back real indoor temp
   T_flow_plot = [T_flow_plot y(2)];
   Q_out = [Q_out Q_HP_new];
   %break;
end
%%
%plotting
% figure;
% hold on;
% plot(1:size(troub_t_in,2),troub_t_in);
% plot(1:size(troub_t_flow,2),troub_t_flow);
% plot(1:size(troub_t_out,2),troub_t_out);
% plot(1:size(troub_q_out,2),troub_q_out);
% legend("tin","tflow","tout","qout")

linesize = 1.2;
T_in_plot_without = T_in_plot;

load("t_in_plot_with.mat")

figure;
subplot(2,1,1);
set(gca, 'FontSize', 12, 'LineWidth', linesize); 
hold on;
grid on;
%correct_x_gaps = linspace(0,sim_days*24,n);
plot(1:size(T_in_plot,2),T_in_plot,'LineWidth',linesize,'Color','g');
plot(1:size(T_in_plot_without,2),T_in_plot_without,'LineWidth',linesize);
ylim([18 24])
yline(21,'--','LineWidth',1.2);
legend("Indoor Temperature - Cost Function with Price Coefficient","Indoor Temperature - Cost Function without Price Coefficient","Target Indoor Temperature",'FontSize', 12, 'Location', 'best')

% making the temp and input datasets interpolatable
xlabel("Data Index (half-hourly intervals)");
ylabel("Temperature (°C)")

subplot(2,1,2);
set(gca, 'FontSize', 12, 'LineWidth', linesize); 
hold on;
grid on;
stairs(1:size(p_buy_price(1+day_of_year:n+day_of_year),2),p_buy_price(1+day_of_year:n+day_of_year),'LineWidth',linesize)
xlabel("Data Index (half-hourly intervals)");
ylabel("Price (p/kWh)");
legend("Electricity Buy Price", 'FontSize', 12, 'Location', 'best')



toc;
%% functions


function [c, ceq] = psi_to_q(d)

    global size_decision; 
    global N;
    global x_size;
    global p;

   

   % range is size-N+1 to end for psi terms

   % forming the three vectors for creating equality constraints, which
   % give the right indexes in the decision variables

   psi_vec = size_decision-N+1:size_decision; % psi is the compressor speed AKA multiplier of COP
   
   x_range = 1:x_size*(N+1);
   indoorT_vec = x_range(1:2:end-1);
   flowT_vec = x_range(2:2:end);

   u_range = x_size*(N+1)+1:x_size*(2*N+1); % picks the u range
   outdoorT_vec = u_range(1:2:end-1); 
   Q_vec = u_range(2:2:end); 

   %setting constraints
   for i = 1:N
        % Q(i) = COP(i)*psi(i)
        ceq(i) = d(Q_vec(i))-d(psi_vec(i))*polyvaln(p,[d(flowT_vec(i)),d(outdoorT_vec(i))]); % this means this term must equal zero, POLYVAL GIVES THE COP
   end
        
    c = [];
end
% ISSUE is with the initial terms and the cost function
function objective = cost_function(d)

    global p_buy_price;
    global size_decision;
    global x_size;
    global N;
    global target_indoor_temp;
    global i;

    x_range = 1:x_size*(N+1);
    psi_vec = size_decision-N+1:size_decision;

   f = [];
        for j = 0:N-1
            new = [p_buy_price(i+j)];
            f = [f new]; 
        end

   objective_compressor = 0.1*sum(d(psi_vec(1:1:N)));
   
   objective = sum((d(1:x_size:N*2-1)-target_indoor_temp).^2)+1*d(psi_vec(1));

end

