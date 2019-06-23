global A B C U T_deep T dx dt flag omega

%% User select noise flag:
% noise_flag = 'deterministic';
% noise_flag = 'additive';
% noise_flag = 'multiplicative';       % defunct: do not use (little impact)
% noise_flag = 'multiplicative_kappa'; % multiply kappa by random number
% noise_flag = 'additive_multiplicative'; % multiply kappa by random number also additive



%% Choice of flag
flag = 'none';
% flag = 'cycle_U';
% flag = 'cycle_T';
% flag = 'Jac';


%% parameters for timestepping

transi = 400*12*30*24*60*60;  % in seconds
% ttot   = 1000*12*30;  % run length in days 
ttot = 30*12*30*24*60*60;
h     = 12*60*60;          % in seconds

%% initial conditions
u_init = -0.5;
T_w_init = 18;
T_e_init = 12;

% parameters from Vallis
A = 1/12/30/24/60/60;
B = 2;
C = 1/4/30/24/60/60;
U = -0.45;
T_deep = 6.65;
T = 12;
dx = 7500*1000;
omega = 2*pi/(12*30*24*60*60);

% h_p   = ( b_p/(kappa*a_p))*(a_p-1);
% h_m   = (-b_m/(kappa*a_m))*(a_m-1);


% parameters from Stone et al
% delta_1 = 1.15*30;    %delay time in days
% delta_2 = 5*delta_1;  %delay time in days
% a       = 1/180;
% b       = 1/120;
% c       = 0.9/138;
% f_o     = 1/365;      % frequency seasonal cycle   day^(-1)
% 
% D       = 0.165;        % magnitude stochastic noise


% noise parameters
% sd = D;
% phi = 0.3;
% e_past = 0;
% s1 = RandStream('mt19937ar','seed',1);
% 
% 
% disp(noise_flag)

%% derived parameters
no_steps_transi = transi/h;
no_steps = ttot/h;

%% initialise variables transient
y_vec_transi = zeros(3,no_steps_transi+no_steps);
t_vec_transi = zeros(1,no_steps_transi+no_steps);
% y_vec = zeros(3,no_steps);
% t_vec = zeros(1,no_steps);
y_init = [u_init;T_w_init;T_e_init];
y_vec_transi(:,1) = y_init;


%% integrate transient and real in one step
% for i = 2:no_steps_transi
%     y_vec_transi(:,i) = RK4(t_vec_transi(i-1),y_vec_transi(:,i-1));
%     t_vec_transi(i) = dt*(i-1);
% %     if mod(j*dt,10*30*12)==0
% %         disp(j*dt/(12*30*24*60*60))
% %     end
% end
% 
% y_vec(:,1) = y_vec_transi(:,end);
% 
% for j = 2:no_steps
%     [y_vec(:,j)] = RK4(t_vec(j-1),y_vec(:,j-1));
%     t_vec(j) = dt*(j-1);
% %     if mod(j*dt,10*30*12)==0
% %         disp(j*dt/(12*30*24*60*60))
% %     end
% end

f1 = @(t,x) [B*(x(3)-x(2))/2/dx-C*(x(1)-U);...
           x(1)*(T_deep-x(3))/2/dx-A*(x(2)-T);...
           x(1)*(x(2)-T_deep)/2/dx-A*(x(3)-T)];
f2 = @(t,x) [B*(x(2)-x(1))/2/dx-C*(x(1)-(1+3*sin(omega*t)));... % Seasonal cycle u*
           x(1)*(T_deep-x(3))/2/dx-A*(x(2)-T);...
           x(1)*(x(2)-T_deep)/2/dx-A*(x(3)-T)];
f3 = @(t,x) [B*(x(3)-x(2))/2/dx-C*(x(1)-U);...
           x(1)*(T_deep-x(3))/2/dx-A*(x(2)-3*sin(omega*t));...
           x(1)*(x(2)-T_deep)/2/dx-A*(x(3)-3*sin(omega*t))];
f4 = @Vallis_ext;
if strcmp(flag,'none')
    for i = 2:no_steps_transi+no_steps
        k1 = f1(t_vec_transi(i-1),y_vec_transi(:,i-1));
        k2 = f1(t_vec_transi(i-1)+h/2,y_vec_transi(:,i-1)+h*k1/2);
        k3 = f1(t_vec_transi(i-1)+h/2,y_vec_transi(:,i-1)+h*k2/2);
        k4 = f1(t_vec_transi(i-1)+h,y_vec_transi(:,i-1)+h*k3);
        y_vec_transi(:,i) = y_vec_transi(:,i-1)+h*(k1+2*k2+2*k3+k4)/6;
        t_vec_transi(i) = h*(i-1);
    end
end

if strcmp(flag,'cycle_U')
    for i = 2:no_steps_transi+no_steps
        k1 = f2(t_vec_transi(i-1),y_vec_transi(:,i-1));
        k2 = f2(t_vec_transi(i-1)+h/2,y_vec_transi(:,i-1)+h*k1/2);
        k3 = f2(t_vec_transi(i-1)+h/2,y_vec_transi(:,i-1)+h*k2/2);
        k4 = f2(t_vec_transi(i-1)+h,y_vec_transi(:,i-1)+h*k3);
        y_vec_transi(:,i) = y_vec_transi(:,i-1)+h*(k1+2*k2+2*k3+k4)/6;
        t_vec_transi(i) = h*(i-1);
    end
end

if strcmp(flag,'cycle_T')
    for i = 2:no_steps_transi+no_steps
        k1 = f3(t_vec_transi(i-1),y_vec_transi(:,i-1));
        k2 = f3(t_vec_transi(i-1)+h/2,y_vec_transi(:,i-1)+h*k1/2);
        k3 = f3(t_vec_transi(i-1)+h/2,y_vec_transi(:,i-1)+h*k2/2);
        k4 = f3(t_vec_transi(i-1)+h,y_vec_transi(:,i-1)+h*k3);
        y_vec_transi(:,i) = y_vec_transi(:,i-1)+h*(k1+2*k2+2*k3+k4)/6;
        t_vec_transi(i) = h*(i-1);
    end
end

if strcmp(flag,'Jac')
    y_vec_transi = zeros(12,no_steps_transi+no_steps);
    y_vec_transi(1:3,1) = y_init;
    for i=1:3 
        y_vec_transi(4*i,1)=1; 
    end
    for i = 2:no_steps_transi+no_steps
        k1 = f4(t_vec_transi(i-1),y_vec_transi(:,i-1));
        k2 = f4(t_vec_transi(i-1)+h/2,y_vec_transi(:,i-1)+h*k1/2);
        k3 = f4(t_vec_transi(i-1)+h/2,y_vec_transi(:,i-1)+h*k2/2);
        k4 = f4(t_vec_transi(i-1)+h,y_vec_transi(:,i-1)+h*k3);
        y_vec_transi(:,i) = y_vec_transi(:,i-1)+h*(k1+2*k2+2*k3+k4)/6;
        t_vec_transi(i) = h*(i-1);
    end
end

y_vec = y_vec_transi(:,no_steps_transi+1:end);
t_vec = t_vec_transi(no_steps_transi+1:end);
figure(1)
plot(t_vec/60/60/24/30/12,y_vec(1,:))
figure(2)
plot(y_vec(1,:),y_vec(3,:)-y_vec(2,:))
% figure(3)
% plot3(y_vec(1,:),y_vec(2,:),y_vec(3,:))
% figure(4)
% plot(t_vec/60/60/24/30/12,y_vec(3,:)-y_vec(2,:))
% figure(5)
% plot(t_vec/60/60/24/30/12,y_vec(2,:))




       