clear; close all;
yalmip('clear')

%% Parameters
t_sim=24;
T=1;
v0=18;
x0=[15;0];


z=zeros(2,t_sim);
z(:,1)=x0;
u=zeros(1,t_sim);

C = 1;
f = 1;
Ki=1;

Anu = exp(-Ki/C*T);
Bu = 0.5;
Bv = (1-exp(-Ki/C*T));
W = 0.1;
F = 0.5;
K=1;
horizon = 10;
A = [Anu,Bv*F;0,F];
B = [Bu; 0];
H = [Bv; 0];
E = [Bv*K; K];
xmax = 22;
xmin = 18;
umax=2;
c=ones(1,horizon);
kappa = 1000;

pred = v0*ones(t_sim/T+horizon,1);
v_real = pred + generateUWhiteNoise([-1 1],t_sim/T+horizon)';
vtil=ones(t_sim/T,1);



%% Variables

M = sdpvar(horizon,horizon);
h = sdpvar(horizon,1);
pred_use = sdpvar(horizon,1);
tstart = sdpvar(1);
xstart = sdpvar(2,1);
slack = sdpvar(horizon,1);

%% Formulation
constraints = [];
for k=1:horizon
    for j= 1:horizon
        constraints = [constraints, M(k,j)*(k>=j)==0];
    end
end
for t=1:horizon
L = [1,0]*A^t;
G = zeros(2,t);
V = zeros(2,t);
W = zeros(2,t);
G(:,t)=B;
V(:,t)=H;
W(:,t)=E;
for k = 1:t-1
    G(:,t-k)=A*G(:,t-k+1);
    V(:,t-k)=A*V(:,t-k+1);
    W(:,t-k)=A*W(:,t-k+1);
end
G = [1,0]*G;
V = [1,0]*V;
W = [1,0]*W;
y = sdpvar(1,t);
for k = 1:t
    co = 0;
    for j = 1:k-1
        co = co + G(k)*M(k,j);
    end
    co = co+ W(k);
    constraints = [constraints, y(k)==co];
end
constraints = [constraints, 1.4*sqrt(y*y')+G*h(1:t)<=xmax-L*xstart-V*pred_use(1:t)+slack(t)];
constraints = [constraints, 1.4*sqrt(y*y')-G*h(1:t)<=-xmin+L*xstart+V*pred_use(1:t)-slack(t)];
constraints = [constraints, slack(t)>=0];
constraints = [constraints, h(t)>=0, h(t)<=umax];
end
cost = c*h+kappa*ones(1,horizon)*slack;

%% Time-dependent implementation
xold=x0;
for k=1:t_sim/T
newconstraints = [constraints, xstart==xold];
newconstraints = [newconstraints, pred_use==pred(k+1:k+horizon)];
assign(M,ones(horizon));
diagnostic = optimize(newconstraints,cost,sdpsettings('solver','fmincon','usex0',1));
u(k)=value(h(1))
% Simulation
z(:,k+1) = A*z(:,k)+B*u(k)+H*(v_real(k));
z(2,k+1) = v_real(k)-pred(k);
xold=z(:,k+1);
end
%% Results
figure
plot(z(1,:))
