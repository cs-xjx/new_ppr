function p = iter_compute(c, Q, p0, kmax)
%ppr_iter get PPR scores by iteration
%   Q: transition matrix
%   kmax: # if the maximum iterations

p = p0;
for iter = 1 : kmax
    p = c * Q * p + p0;
end
% p = (1 - c) * p;

end

