%% Set up stochastic noise
phi = 0.99;
e_past = 0;
s1 = RandStream('mt19937ar','seed',1);
x = zeros(1,10000);

%%
for i = 1:1:10000
    %stochastic term
    e_stoch = phi*e_past + ((1-phi^2)^0.5)*randn(s1,1,1);
    e_past = e_stoch;
    x(i) = e_stoch;
end

plot(x)