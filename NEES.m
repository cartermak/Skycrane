function [eps_xk] = NEES(x_err, P_kp)

eps_xk = x_err'* P_kp^-1 *x_err;


end