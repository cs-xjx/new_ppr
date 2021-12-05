function [res] = BearS_uncertain_mut(a, c, src, tar, qu, ds)
time_total = tic;
fprintf('\n >> Start bearS_MutualExclusion ... \n');

addPathes;
debug = false;

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
a = flatten_mut(a, src, tar0, n_certain);
% Q = flattenUnGraph_mut(A, src, tar, na);
period_flatten = toc(time_flatten);

% pre-process
%
% 2. pre-process
%
time_preprocess = tic;
fprintf(' >> 2. Pre-process ... \n');
[PRE_RET] = BearS_Pre(a, 1 - c, 0, debug);
period_preprocess = toc(time_preprocess);


%
% 3. online query 
%
time_online_query = tic;  
fprintf(' >> 3. Online Query ... \n');
res = cell(qu_num, 1);
for qi = 1 : qu_num
    fprintf(' query #%d\n', qi);
    res{qi} = BearS_Query(qu{qi}, 1 - c, PRE_RET);
%     [x, y] = sort(res{qi}, 'descend');
%     x(1 : 10)
%     y(1 : 10)
end
period_online_query = toc(time_online_query);

period_total = toc(time_total);

fprintf('\n========== BearS_uncertain_mut ============\n');

fprintf(' >> Dataset                                 :  %s \n',   ds);
fprintf(' >      # of nodes (certain)                :  %d \n',   n_certain);
fprintf(' >      # of edges (certain)                :  %d \n',   m_certain);
% fprintf(' >      # of possible worlds                :  %d \n',   npw);
fprintf(' >      # of queries                        :  %d \n\n', qu_num);

fprintf(' >> Parameters\n');
fprintf(' >      decay factor (c)                    :  %d \n',   c);
% fprintf(' >      low-rank SVD (r)                    :  %d \n',   r);
% fprintf(' >      # of partitions                     :  %d \n\n', n_blocks);

fprintf(' >> Total Time                              :  %fs \n', period_total);
fprintf(' >    1) Preprocessing                      :  %f s \n', period_flatten + period_preprocess);
fprintf(' >         a1) Flatten uncertain graph      :  %f s \n', period_flatten);
fprintf(' >         a2) pre-process                  :  %f s \n', period_preprocess);
fprintf(' >    2) Online Query (# of queries = %d)   :  %f s \n\n', qu_num, period_online_query);

fprintf(' >> Average Amortised Total Time per Query              :  %fs \n', period_total / qu_num);
fprintf(' >> Average Amortised Preprocessing Time per Query      :  %fs \n', (period_flatten + period_preprocess) / qu_num);

whos
