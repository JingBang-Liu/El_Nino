function out = hdot(t_now,h_now,h_past_1, h_past_2,e_stoch)
% function to calculate gradient of DO equation

global a b c f_o D_add D_mult noise_flag kappa kappa_stoch 

% out = (1+e_stoch)*T_now-T_now^3-alpha*T_delay;

if strcmp(noise_flag,'additive')
    kappa_stoch = kappa;
    out = a*A_h(h_past_1)-b*A_h(h_past_2) + ...
              c*cos(2*pi*f_o*t_now)+D_add*c*e_stoch;
elseif strcmp(noise_flag,'multiplicative')
    kappa_stoch = kappa;
    out = (1+D_mult*c*e_stoch)*(a*A_h(h_past_1)-b*A_h(h_past_2)) + ...
              c*cos(2*pi*f_o*t_now);
elseif strcmp(noise_flag,'multiplicative_kappa')
    kappa_stoch = kappa*(1+D_mult*c*e_stoch);
    out = a*A_h(h_past_1)-b*A_h(h_past_2) + ...
              c*cos(2*pi*f_o*t_now);
elseif strcmp(noise_flag,'deterministic')
    kappa_stoch = kappa;
    out = a*A_h(h_past_1)-b*A_h(h_past_2) + ...
              c*cos(2*pi*f_o*t_now);
elseif strcmp(noise_flag,'additive_multiplicative')
    kappa_stoch = kappa*(1+D_mult*c*e_stoch);
    out = a*A_h(h_past_1)-b*A_h(h_past_2) + ...
              c*cos(2*pi*f_o*t_now)+D_add*c*e_stoch; 
else
    error('wrong noise flag')
end

end


