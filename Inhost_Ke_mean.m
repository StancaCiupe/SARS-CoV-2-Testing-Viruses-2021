% Solving Inhost model- based on Ke et al. Mean in table S6
clc
clear all

global beta1 c k1 beta2 k2 pi1 pi2 Gamma delta10 delta20 w


% %mean values Ke Bioarxiv
c=10;
k1=4;
k2=4;
Gamma=0.01;
beta1=51.35e-8;
delta10=1.98;
pi1=49.8;
beta2=70.67e-8;
delta20=0.53;
pi2=0.34;

% Patient Data

T10 = 4.8e+6;
E10=0;
I10 = 10;
V10=0;
T20=4.8e+8;
E20=0;
I20=1;
V20=0;
y0 = [T10 E10 I10 V10 T20 E20 I20 V20];
Tf = 20;

options2=odeset('NonNegative',1);

[t, y] = ode15s('covid_two_Ke', [0 Tf], y0);
V1 = log10(y(:,4));
V2 = log10(y(:,8));
figure(1)
plot(t,V1, 'r')
hold on
plot(t,V2, 'b')
axis([0 20 0 8])

for k=1:length(y(:,4))
        if round(log10(y(k,4)))==2
            disp(t)
        end
end