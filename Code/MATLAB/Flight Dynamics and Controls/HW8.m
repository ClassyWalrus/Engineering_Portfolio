
[t,x] = ode45(@(t,x) ODEFUN(t,x),[0,4],[2,1]);

X_1 = x(49,1)
X_2 = x(49,2)

function XStep = ODEFUN(t,x)
    XStep = zeros(2,1);
    XStep(1) = -0.6*x(1)-0.7*x(2)+1;
    XStep(2) = -0.7*x(1)-1.65*x(2)+0;
end