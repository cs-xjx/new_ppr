function [res] = flatPPR_MultipleEdge(a, c, kmax, src, tar, qu, ds) 
time_total = tic;
fprintf('\n >> Start flatPPR_MultipleEdge ... \n');

n_certain = size(a, 1); % n_certain: # of nodes in graph
m_certain = nnz(a);

qu_num = numel(qu);
% nsrc = numel(src);

n = max([n_certain, tar{:}]);
tar0 = cellfun(@(x) [0 x], tar, 'UniformOutput', false);

a(n, n) = 0;

%
% 1. flatten uncertain graph
%
time_flatten = tic;
fprintf(' >> 1. Flatten Uncertain Graph ... \n');
% flatten each uncertain edge into certain edges
% Q2 = flattenUnGraph_mul(a, src, tar, n_certain);
a2 = flatten_mul(a, src, tar0, n_certain);
sumq = sum(a2, 2);
d = 1 ./ sumq;
d(~isfinite(d)) = 0;
Q = a2' * spdiags(d, 0, n, n);
clear a2;
period_flatten = toc(time_flatten);

%
% 2. online query 
%
time_online_query = tic;  
fprintf(' >> 2. Online Query ... \n');
res = cell(qu_num, 1);
for qi = 1 : qu_num
    fprintf(' query #%d\n', qi);
    res{qi} = ppr_iter(c, Q, qu{qi}, kmax);
end
period_online_query = toc(time_online_query);

period_total = toc(time_total);

fprintf('\n========== flatPPR_MultipleEdge ============\n');

fprintf(' >> Dataset                                 :  %s \n',   ds);
fprintf(' >      # of nodes (certain)                :  %d \n',   n_certain);
fprintf(' >      # of edges (certain)                :  %d \n',   m_certain);
% fprintf(' >      # of possible worlds                :  %d \n',   npw);
fprintf(' >      # of queries                        :  %d \n\n', qu_num);

fprintf(' >> Parameters\n');
fprintf(' >      decay factor (c)                    :  %d \n',   c);
fprintf(' >      # of iterations                     :  %d \n', kmax);
% fprintf(' >      low-rank SVD (r)                    :  %d \n',   r);
% fprintf(' >      # of partitions                     :  %d \n\n', n_blocks);

fprintf(' >> Total Time                              :  %fs \n', period_total);
fprintf(' >    1) Preprocessing - flatten            :  %fs \n', period_flatten);
fprintf(' >    2) Online Query (# of queries = %d)   :  %fs \n\n', qu_num, period_online_query);

fprintf(' >> Average Amortised Total Time per Query              :  %fs \n', period_total / qu_num);
% fprintf(' >> Average Amortised Total Time Time per Query per PW  :  %fs \n', period_total / qu_num / npw);
fprintf(' >> Average Amortised Preprocessing Time per Query      :  %fs \n', period_flatten / qu_num);

whos
