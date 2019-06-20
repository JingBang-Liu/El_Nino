

function [y_n] = RK4(t,y)
global A B C U T_deep T dx dt flag omega
f1 = @(t,x) [B*(x(3)-x(2))/dx-C*(x(1)-U);...
           x(1)*(T_deep-x(3))/2/dx-A*(x(2)-T);...
           x(1)*(x(2)-T_deep)/2/dx-A*(x(3)-T)];
f2 = @(t,x) [B*(x(2)-x(1))/dx-C*(x(1)+(1+2*sin(omega*t)));... % Seasonal cycle u*
           x(1)*(T_deep-x(3))/2/dx-A*(x(2)-T);...
           x(1)*(x(2)-T_deep)/2/dx-A*(x(3)-T)];
       
% if strcmp(flag,'none')
%     k1 = dt*f1(t,y);
%     k2 = dt*f1(t,y+k1./2);
%     k3 = dt*f1(t,y+k2./2);
%     k4 = dt*f1(t,y+k3);
% end
% if strcmp(flag,'anual_cycle')
%     k1 = dt*f2(t,y);
%     k2 = dt*f2(t+dt/2,y+k1./2);
%     k3 = dt*f2(t+dt/2,y+k2./2);
%     k4 = dt*f2(t+dt,y+k3);
% end
% y_n = y+(k1+2*k2+2*k3+k4)/6;

if strcmp(flag,'none')
    k1 = f1(t,y);
    k2 = f1(t+0.5*dt,y+dt*0.5*k1);
    k3 = f1(t+3*dt/8,y+dt*3*k1/32+dt*9*k2/32);
    k4 = f1(t+12*dt/13,y+dt*1932*k1/2197-dt*7200*k2/2197+dt*7296*k3/2197);
    k5 = f1(t+dt,y+dt*439*k1/216-dt*8*k2+dt*3680*k3/513-dt*845*k4/4104);
end
if strcmp(flag,'anual_cycle')
    k1 = dt*f2(t,y);
    k2 = dt*f2(t+0.5*dt,y+0.5*k1);
    k3 = dt*f2(t+3*dt/8,y+3*k1/32+9*k2/32);
    k4 = dt*f2(t+12*dt/13,y+1932*k1/2197-7200*k2/2197+7296*k3/2197);
    k5 = dt*f2(t+dt,y+439*k1/216-8*k2+3680*k3/513-845*k4/4104);
end
y_n = y+dt*(25*k1/216+1408*k3/2565+2197*k4/4104-k5/5);
end
