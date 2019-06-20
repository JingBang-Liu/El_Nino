function out = A_h(h_in)

global a_p a_m b_p b_m kappa_stoch %h_p h_m

h_p_tmp   = ( b_p/(kappa_stoch*a_p))*(a_p-1);
h_m_tmp   = (-b_m/(kappa_stoch*a_m))*(a_m-1);


if h_in > h_p_tmp
    
   out = b_p + (b_p/a_p)*(tanh((kappa_stoch*a_p/b_p)*(h_in-h_p_tmp))-1);
    
elseif h_in >= h_m_tmp
    
   out = kappa_stoch*h_in;
    
elseif h_in < h_m_tmp
    
   out = -b_m - (b_m/a_m)*(tanh((kappa_stoch*a_m/b_m)*(h_m_tmp-h_in))-1);
    
else
    error('wrong!')
end

end