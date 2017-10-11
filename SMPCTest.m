v0=18;
t_sim=24;
T=1;
horizon=10;
pred = v0*ones(t_sim/T+horizon,1);
for t=1:t_sim/T+horizon
    pred(t)=v0-2*sin(2*pi*t*T/t_sim);
end
v_real = pred + generateUWhiteNoise([-1 1],t_sim/T+horizon)';
[SMPC18,uS,ruleS]=SMPCTrue(pred,v_real,v0);
[DMPC18,uD,ruleD]=DMPC(pred,v_real,v0);
[PB18,uP,ruleP]=PB(v_real,v0);

figure
hold on
plot(SMPC18,'b')
plot(DMPC18,'r')
plot(PB18,'g')
plot([0,25],[18,18],'k')
plot([0,25],[22,22],'k')
legend('SMPC','DMPC','PB','Comfort bound')
xlabel('Time [h]')
ylabel('Temperature [°C]')

figure
hold on
plot(SMPC18,'b')
plot(pred(1:t_sim/T),'c')
plot(v_real(1:t_sim/T),'m')
plot([0,25],[18,18],'k')
plot([0,25],[22,22],'k')
legend('SMPC','Weather prediction','True outer temperature','Comfort bound')
xlabel('Time [h]')
ylabel('Temperature [°C]')