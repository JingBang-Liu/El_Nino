
steps = 10001;
TT = linspace(-5,10,steps);
TT_max = 0;
UU_max = 0;

%% Choice of which cycle to use
flag = 'none';
% flag = 'cycle_U';
% flag = 'cycle_T';

%% Set the length of time interval
transi = 8*12*30*24*60*60;  % in seconds, this is 8 years
ttot = 200*12*30*24*60*60; % this is 200 years
h     = 12*60*60;          % in seconds, this is 0.5 day

%% initial conditions
u_init = -0.5; % in m*sec^(-1)
T_w_init = 18; % in C
T_e_init = 12; % in C

for j = 1:steps
    % parameters from Vallis
    A = 1/12/30/24/60/60; % sec^(-1)
    B = 2; % m^2 * sec^(-2) * C^(-1)
    C = 1/4/30/24/60/60; % sec^(-1)
    U = -0.45; % this is u* in the paper, in m*sec^(-1)
    T_deep = TT(j); % this is T bar in the paper, in C
    T = 12; % this is T* in the paper, in C. This is different from the paper, as in the paper Vallis used 12.
    dx = 7500*1000; % m
    omega = 2*pi/(12*30*24*60*60); % the parameter for annual cycle, as in seconds a year

    %% derived parameters
    no_steps_transi = transi/h; % total timesteps: transient timesteps and timesteps
    no_steps = ttot/h; % timesteps 

    %% initialise variables transient
    y_vec_transi = zeros(3,no_steps_transi+no_steps);
    t_vec_transi = zeros(1,no_steps_transi+no_steps);
    y_vec = zeros(3,no_steps);
    t_vec = zeros(1,no_steps);
    y_init = [u_init;T_w_init;T_e_init];
    y_vec_transi(:,1) = y_init;


    %% Using RK4 


    f1 = @(t,x) [B*(x(3)-x(2))/2/dx-C*(x(1)-U);...  % ODE function without seasonal cycle
               x(1)*(T_deep-x(3))/2/dx-A*(x(2)-T);...
               x(1)*(x(2)-T_deep)/2/dx-A*(x(3)-T)];
    f2 = @(t,x) [B*(x(2)-x(1))/2/dx-C*(x(1)-(1+3*sin(omega*t)));... % Seasonal cycle on u*
               x(1)*(T_deep-x(3))/2/dx-A*(x(2)-T);...
               x(1)*(x(2)-T_deep)/2/dx-A*(x(3)-T)];
    f3 = @(t,x) [B*(x(3)-x(2))/2/dx-C*(x(1)-U);...
               x(1)*(T_deep-x(3))/2/dx-A*(x(2)-3*sin(omega*t));...  % Seasonal cycle on T*
               x(1)*(x(2)-T_deep)/2/dx-A*(x(3)-3*sin(omega*t))];


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

    y_vec = y_vec_transi(:,no_steps_transi+1:end);
    t_vec = t_vec_transi(no_steps_transi+1:end);

    %% Bifurcation diagram
    u = y_vec(1,:);
    z = z1test(u(100000:200:end));
    if z>0.1
        
    gap = 2;
    k = 1;
    local_max = 0;
    for i = 1+gap:gap:no_steps-gap
        if u(i)>u(i-gap)&&u(i)>u(i+gap)
            local_max(k) = u(i);
            k = k+1;
        end
    end
    s = size(local_max);
    s = s(2);
    TT_temp = ones(1,s)*TT(j);
    UU_max = [UU_max,local_max];
    TT_max = [TT_max,TT_temp];
    j
    end
end

scatter(TT_max,UU_max)