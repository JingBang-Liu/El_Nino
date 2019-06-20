% RK4.m

function [tout,hout,enoise] = RK4(t_now,hvec)

global dt delta_1 delta_2 phi e_past s1

t_index     = round(t_now/dt);
delta_index_1 = round(delta_1/dt);
delta_index_2 = round(delta_2/dt);

%stochastic term
e_stoch = phi*e_past + ((1-phi^2)^0.5)*randn(s1,1,1);
e_past = e_stoch;

% e_stoch = phi*e

hnow   = hvec(t_index);
hdelay_1 = hvec(t_index-delta_index_1);
hdelay_2 = hvec(t_index-delta_index_2);

k1 = hdot(t_now     ,hnow            ,hdelay_1,hdelay_2,e_stoch);
k2 = hdot(t_now+dt/2,hnow + (dt/2)*k1,hdelay_1,hdelay_2,e_stoch);
k3 = hdot(t_now+dt/2,hnow + (dt/2)*k2,hdelay_1,hdelay_2,e_stoch);
k4 = hdot(t_now+dt  ,hnow + dt*k3    ,hdelay_1,hdelay_2,e_stoch);

tout = t_now+dt;
hout = hnow + (dt/6)*(k1+2*k2+2*k3+k4);
enoise = e_stoch;