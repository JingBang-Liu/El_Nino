A = 1/12/30/24/60/60;
B = 2;
C = 1/4/30/24/60/60;
U = -0.45;
T_avg = 0;
T = 3;
dx = 7500*1000;
omega = 2*pi/(12*30*24*60*60);
x_init = [-0.5;15;10];
f1 = @(t,x) [B*(x(2)-x(1))/2/dx-C*(x(1)-U);... % Normal u*
           x(1)*(T_avg-x(3))/2/dx-A*(x(2)-T);...
           x(1)*(x(2)-T_avg)/2/dx-A*(x(3)-T)];
       
f2 = @(t,x) [B*(x(2)-x(1))/dx-C*(x(1)-(1+3*sin(omega*t)));... % Seasonal cycle u*
           x(1)*(T_avg-x(3))/2/dx-A*(x(2)-T);...
           x(1)*(x(2)-T_avg)/2/dx-A*(x(3)-T)];
f3 = @(t,x) [B*(x(3)-x(2))/dx-C*(x(1)-U);...
           x(1)*(T_deep-x(3))/2/dx-A*(x(2)-3*sin(omega*t));...
           x(1)*(x(2)-T_deep)/2/dx-A*(x(3)-3*sin(omega*t))];
       
[t,a] = ode45(f1,[0:60*60*12:48*12*30*24*60*60],x_init);
figure(1)
plot(t(4000:1:end)/(60*60*24*30*12),a(4000:1:end,1))
figure(2)
plot(a(3000:1:end,1),a(3000:1:end,3)-a(3000:1:end,2))
% figure(3)
% plot(t(4000:1:end)/(60*60*24*30*12),a(4000:1:end,3)-a(4000:1:end,2))
% plot(t(193:end),a(3*12*30*24*60*60:1:end,1))
% figure(2)
% plot(a(:,1),a(:,3)-a(:,2))