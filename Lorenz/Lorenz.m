h = 0.001;
T = 100;
iterations = T/h;
x_init = [1;1;1];
x_vec = zeros(3,iterations+1);
x_vec(:,1) = x_init;
t_vec = linspace(0,T,iterations+1);
sigma = 10;
rho = 28;
beta = 8/3;

%% Set up stochastic noise
phi = 0;
e_past = 0;
s1 = RandStream('mt19937ar','seed',1);
D = 0;

%% The function we want to solve
f = @(t,x,stoch) [-sigma.*x(1)+sigma.*x(2)+D*stoch;...
                 x(1).*(rho-x(3))-x(2)+D*stoch;...
                 x(1).*x(2)-beta.*x(3)+D*stoch];
        
        
for i = 2:1:iterations+1
    %stochastic term
    e_stoch = phi*e_past + ((1-phi^2)^0.5)*randn(s1,1,1);
    e_past = e_stoch;
    k1 = f(t_vec(i-1),x_vec(:,i-1),e_stoch);
    k2 = f(t_vec(i-1)+h/2,x_vec(:,i-1)+h*k1/2,e_stoch);
    k3 = f(t_vec(i-1)+h/2,x_vec(:,i-1)+h*k2/2,e_stoch);
    k4 = f(t_vec(i-1)+h,x_vec(:,i-1)+h*k3,e_stoch);
    x_vec(:,i) = x_vec(:,i-1) + h*(k1+2*k2+2*k3+k4)/6;
end

plot3(x_vec(1,:),x_vec(2,:),x_vec(3,:))


% function [tt,y] = RK4(t,x,h)
%     sigma = 10;
%     rho = 28;
%     beta = 8/3;
%     f = @(t,x) [-sigma.*x(1)+sigma.*x(2);...
%                 x(1).*(rho-x(3))-x(2);...
%                 x(1).*x(2)-beta.*x(3)];
%     k1 = f(t,x);
%     k2 = f(t+h/2,x+h*k1/2);
%     k3 = f(t+h/2,x+h*k2/2);
%     k4 = f(t+h,x+h*k3);
%     y = x + h*(k1+2*k2+2*k3+k4)/6;
%     tt = t + h;
% end