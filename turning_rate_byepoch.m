function [alpha, pe] =turning_rate_byepoch(r,X)
% construction of series of turning rates alpha_1, ... alpha_[T/dim epoch]
% epoch = r
T=length(X);

for t=1:floor(T/r)
    Y=X((t-1)*r+1 :t*r);
    bit_ord=lag_patterns_length3(Y,1);
    alpha(t)=turning_rate(bit_ord);
   % pe(t)=PE(bit_ord);
end
pe=1;
end 