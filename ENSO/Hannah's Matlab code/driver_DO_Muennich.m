% driver_DO_Muennich.m

clear all

% deterministic delayed oscillator model for ENSO

% following Muennich, Cane and Zebiak 1991 JAS
% and Stone, Saparin, Huppert and Price, 1998, GRL.


global a b c f_o D_add D_mult a_p a_m b_p b_m kappa kappa_stoch h_p ...
      h_m dt delta_1 delta_2 phi e_past s1 noise_flag


%% User select noise flag:
% noise_flag = 'deterministic';
noise_flag = 'additive';
% noise_flag = 'multiplicative';       % defunct: do not use (little impact)
% noise_flag = 'multiplicative_kappa'; % multiply kappa by random number
% noise_flag = 'additive_multiplicative'; % multiply kappa by random number also additive
  
%% variables
% h_vec      thermocline depth as a function of time
% T_now      Temperature now
%
% t_vec      time - continuous - in days
% t_now      time - continuous - in days
% dt         time - step       - in days
% t_tot      time - total      - in days
% 
% t_i        time - index
% no_steps   no. time steps


%% parameters for timestepping

transi = 100*12*30;  % in days
% ttot   = 1000*12*30;  % run length in days 
% ttot   = 500*12*30;  % run length in days %% testing
ttot = 300*12*30;
dt     = 0.5;          % in days
h_init = 0.25;

% parameters from Tziperman et al, 1994, Science
a_p   = 1;    %must be > 1
a_m   = 1;    %must be > 1
b_p   = 1;
b_m   = b_p/5;    %must be less than b_p

kappa         = 1.6;
kappa_stoch   = kappa;

h_p   = ( b_p/(kappa*a_p))*(a_p-1);
h_m   = (-b_m/(kappa*a_m))*(a_m-1);


% parameters from Stone et al
delta_1 = 1.15*30;    %delay time in days
delta_2 = 5*delta_1;  %delay time in days
a       = 1/180;
b       = 1/120;
c       = 0.9/138;
f_o     = 1/365;      % frequency seasonal cycle   day^(-1)

D_add       = 0.165;        % magnitude of additive stochastic noise
D_mult      = 40;          % magnitude of multipicative stochastic noise

% noise parameters
phi = 0;
e_past = 0;
s1 = RandStream('mt19937ar','seed',1);


disp(noise_flag)

%% derived parameters
no_steps = ttot/dt;
no_steps_transi = transi/dt;

%% initialise variables transient
% need to initialise backawrd in time dep on delay - (delta_2 > delta_1)
delta_index = floor(delta_2/dt)+1;
h_vec_transi = zeros(no_steps_transi,1);
t_vec_transi = zeros(no_steps_transi,1);
for i = 1:delta_index
    h_vec_transi(i) = h_init;
end

%% integrate transient and real in one step
% disp('Integrate to remove transient')
for i = delta_index:no_steps_transi+no_steps
    t_now = dt*(i);    
    [t_vec_transi(i+1),h_vec_transi(i+1)] = RK4(t_now,h_vec_transi);
    if mod(i*dt,10*30*12)==0
        disp(i*dt/(30*12))
    end
end

h_vec = h_vec_transi(no_steps_transi+1:end);
t_vec = (0:dt:ttot)';

%% smooth timeseries using 12 months moving average

no_smooth_over = (12*30/dt);  % in days

h_vec_smooth=zeros(length(h_vec)-(no_smooth_over-1),1);
t_vec_smooth=zeros(length(t_vec)-(no_smooth_over-1),1);
for i_sm = 1:no_smooth_over
    h_vec_smooth = h_vec_smooth+h_vec(i_sm:end-(no_smooth_over-i_sm));
    t_vec_smooth = t_vec_smooth+t_vec(i_sm:end-(no_smooth_over-i_sm));
end
h_vec_smooth = h_vec_smooth/no_smooth_over;
t_vec_smooth = t_vec_smooth/no_smooth_over;
h_vec_smooth = h_vec_smooth(no_smooth_over:end-11);
t_vec_smooth = t_vec_smooth(no_smooth_over:end-11);

% 
figure
plot(t_vec/(12*30),h_vec)
xlabel('Time / years')
ylabel('Height / meters')
axis([0 90 -0.6 0.6])
hold on
plot(t_vec_smooth/(12*30),h_vec_smooth,'k')

figure
plot(t_vec(1000:1000:end)/(12*30),h_vec(1000:1000:end))
hold on
plot(t_vec_smooth(1000:1000:end)/(12*30),h_vec_smooth(1000:1000:end))

%% Test for chaos

Z = z1test(h_vec(1:600:end))

z = z1test(h_vec_smooth(1:1000:end))

%% Rcurrence Plot
% crp(h_vec(10000:20000))
% crp(h_vec_smooth(10000:20000))

%% save desired variables
clear('e_past','h_init','h_vec_transi','i','i_sm','no_smooth_over','s1','t_now','t_vec_transi','transi','ttot')

% save(['100415_DOtimeseries_D',num2str(D),'_kappa',num2str(kappa),'_',noise_flag,'.mat'])










