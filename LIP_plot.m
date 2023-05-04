clear variables;
close all;
H=0.53;
Tc=sqrt(H/9.8);

x0=-0.1;
dx0=0.2;
td=0.4;

tLind=linspace(0,td,1000);
xt=x0*cosh(tLind./Tc)+Tc*dx0*sinh(tLind./Tc);
dxt=x0/Tc*sinh(tLind./Tc)+dx0*cosh(tLind./Tc);

figure();
subplot(2,1,1)
plot(tLind,xt);
ylabel('xt');
subplot(2,1,2)
plot(tLind,dxt);
ylabel('dxt');
