function res = rlog(log_vals)
% function RLOG computes log(a0+a1+...+an)
% input form : log_vals = [log(a0),log(a1),...,log(an)]
% assumptions: log(a0) > log(a1) > ... > log(an)
log_vals = sort(log_vals,'descend');

n_inf    = sum(isinf(log_vals));
if n_inf > 1
    log_vals((end-n_inf+1):end)=-10000;
end

if length(log_vals)> 2
    res = log_vals(1)+log(1+exp(rlog(log_vals(2:end))-log_vals(1)));
elseif length(log_vals) == 2
    res = log_vals(1)+log(1+exp(log_vals(2)-log_vals(1)));
else
    res = log_vals;
end