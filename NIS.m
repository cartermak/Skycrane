function [eps_yk] = NIS(y_err, Sk)


eps_yk = y_err' * Sk^-1 * y_err;


end