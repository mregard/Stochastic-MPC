Uout = 0.18;
Uin  = 0.6;

L    = 5;
W    = 4;
Cw   = 2;
Cl   = 4*L;
H    = 2.4;

Tout = 273.15;

T1ref = 273.15+20;
T2ref = 273.15+22;
T3ref = 273.15+24;
T4ref = 273.15+26;

Cin = 1005;

[temp, date, time] = GetData('mybetemp2017.txt');

temp1 = temp(1:(7*24))';
u = 1:length(temp1);

%P = 8.19003354001708;;
%I = 0.0333545530895414;
%D = 78.5648800912167;
P = 30;
I = 1;
D = 0;
Ts = 1/69;
Tsamp = 1;
test1 = temp(1)*ones(1,60);
for x = 2:length(temp1)
    
     test1 = [test1 temp(x)*ones(1,60)];
end

time = 0:1:length(test1)-1;
var.time = time;
var.signals.values = test1'+273.15;
var.signals.dimension = 1;

%%
simPWR = pwrSINGLE.signals.values;

summ = 0;
for i = 1:10081
    summ = summ + simPWR(i);
end
summ = summ/length(simPWR)