clc
close all
clear all
yalmip('clear')

%general
t_sim   = 24;  % simulation time (hours)
T       = 1;    % sampling time (hours)
v0      = 20;  % initial ouside temperature
x0      = 15;  % initial inside temperature
x_des   = 20;   % desired room temperature

%model
C = 10;         % thermal capacity of room
K = 1;          % heat transfer coefficient
F = 0;          % predicted error autoregressive model
k = 1;          % predicted error autoregressive model

%MPC
H       = 10;   %prediction horizon
nx      = 2;    %number of states 
nu      = 1;    %number of inputs
u_min   = -10;   % minimal control input (heat/h)
u_max   = 10;    % max heat per hour
Q       = 1;
R       = 1;

%simulation
w           = randn(t_sim/T,1);  %random disturbance on outside temperature
x           = zeros(t_sim/T,1);         % room temperature
u           = zeros(t_sim/T,1);         % room temperature
v_bar       = ones(t_sim/T + H,1) *20;      % outside temp forecast
v           = zeros(t_sim/T + H,1);         % outside temp real
v_tilde     = zeros(t_sim/T,1);         % predicted error of outside temperature
x(1)        = x0;
v(1)        = v0;
v_tilde(1)  = 0;

%!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
v(2:end) = v0; %assume just constant ouside temp for a start..
%!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

z   =   [x, v_tilde];
z_old = z(1,:);
for i=2:1:t_sim/T
%    v_tilde(i) = v_tilde(i-1) + k*w(i);

    z_mpc = sdpvar(H+1,nx); %state as variable to be optimized (state for all times of horizon H)
    u_mpc = sdpvar(H, nu);
    constraints = [];
    constraints = [constraints, z_mpc(1,:) == z_old];  %x_old is starting state for each time step
    cost =0;
    for j=1:1:H  %MPC
        cost = cost + (x_des - Q*z_mpc(j+1,1))^2 + R*u_mpc(j)^2;
        constraints = [constraints, z_mpc(j+1,1)== z_mpc(j,1) + 1/C*(K*(v_bar(i+j)) - z_mpc(j,1)) + u_mpc(j)];  %MPC knows prediction but not prediction error (v_tilde)
        constraints = [constraints, z_mpc(j+1,2)== z_mpc(j,2)*F];
        constraints = [constraints, u_min <= u_mpc(j) <= u_max];        
    end
    
    options = sdpsettings('verbose', 0 ,'solver', 'quadprog', 'savedebug', 0, 'savesolveroutput',0 );
%     options = sdpsettings('verbose', 0 ,'solver', 'fmincon', 'savedebug', 0, 'savesolveroutput',0, 'usex0', 1 );
    diagnostic = optimize(constraints, cost, options); %RUN MPC potimizer at each time step i

    z_horizon = double(z_mpc); %predicted states by MPC
    u_horizon = double(u_mpc); %predicted (optimized) inputs by MPC
    u(i) = u_horizon(1,:);  %only apply first input
    
    %simulate
    v(i)       = v_bar(i) + w(i);                               %outside temp
    x(i)       = x(i-1) + 1/C*(K*(v(i-1) - x(i-1)) + u(i));     %insie temp
    
    z_old = [x(i), 0]
    
end
%%
figure
tt=1:t_sim/T;
plot(tt,x_des*ones(length(tt),1),tt,x)
legend('temperature ref','temperature actual')
xlabel('time [h]')
ylabel('room temperature [°C]')
ylim([15 25])

