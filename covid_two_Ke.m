function xdot = covid_two_Ke(t,x)

global beta1 c k1 beta2 k2 pi1 pi2 Gamma delta10 delta20

xdot = zeros(7,1);
xdot(1) = -beta1*x(1)*x(4);
xdot(2) = beta1*x(1)*x(4)-k1*x(2);
xdot(3)=k1*x(2)-delta10*x(3);
xdot(4)=pi1*x(3)-c*x(4);
xdot(5) = -beta2*x(5)*x(8);
xdot(6) = beta2*x(5)*x(8)-k2*x(6);
xdot(7)=k2*x(6)-delta20*x(7);
xdot(8)=pi2*x(7)-c*x(8)+Gamma*x(4);
