function [theta_p, theta_n]= getsoc(c_ss_p,c_ss_n,p)

theta_p = c_ss_p/p.c_s_p_max;
theta_n = c_ss_n/p.c_s_n_max;



end