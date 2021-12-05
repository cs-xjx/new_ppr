function p = ppr_iter(c, Q, qu, kmax)
%ppr_iter get PPR scores by iteration
%   Q: transition matrix
%   qu: query nodes
%   kmax: # if the maximum iterations

%get PPR
n = size(Q, 1); % n: # of nodes in graph
s = sparse(qu, 1, 1 / numel(qu), n, 1); % select nodes in qu as important nodes
p = s;
for iter = 1 : kmax
    p = c * Q * p + s;
end
p = (1 - c) * p;

end

