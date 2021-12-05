function [res] = flatApxPPR_MutualExclusion(a, c, n_blocks, src, tar, qu, r, ds)
time_total = tic;
fprintf('\n >> Start flatApxPPR_MutualExclusion ... \n');


n_certain = size(a, 1);
m_certain = nnz(a);

qu_num = numel(qu);
% nsrc = numel(src); % nsrc: # of uncertain edge

n = max([n_certain, tar{:}]);
tar0 = cellfun(@(x) [0 x], tar, 'UniformOutput', false);

a(n, n) = 0;

% inserts the newly added node into the graph
% new_nodes = cellfun(@(x) x(x > na), tar, 'UniformOutput', false);
% new_nodes = cat(2, new_nodes{:});
% new_nodes_num = 0;
% while numel(new_nodes) ~= 0
%     node = min(new_nodes);
%     new_nodes_num = new_nodes_num + 1;
%     for i = 1 : numel(tar)
%         tar{i}(tar{i} == node) = na + new_nodes_num;
%     end
%     new_nodes = new_nodes(new_nodes ~= node);
% end
% if new_nodes_num > 0
%     A = blkdiag(A, zeros(new_nodes_num));
% end
% n = na + new_nodes_num;

%
% 1. flatten uncertain graph
%
time_flatten = tic;
fprintf(' >> 1. Flatten Uncertain Graph ... \n');
% flatten each uncertain edge into certain edges
% Q = flattenUnGraph_mut(A, src, tar, n_certain);
a2 = flatten_mut(a, src, tar0, n_certain);
sumq = sum(a2, 2);
d = 1 ./ sumq;
d(~isfinite(d)) = 0;
Q = a2' * spdiags(d, 0, n, n);
clear a2;
period_flatten = toc(time_flatten);

%
% 2. partition flattened pw(Q)
%
time_partition_svd = tic;
fprintf(' >> 2. Partition Collapsed PW ... \n');
if n_blocks > n
    error('the value of k should not greater than n');
end
% [map, ~] = metismex('PartGraphRecursive', Q, n_blocks); % small n_blocks
[map, ~] = metismex('PartGraphKway', Q, n_blocks);     % large n_blocks
[sorted_map, node_idx] = sort(map);
block_size = diff(find(diff([-1 sorted_map n])));
sorted_Q = Q(node_idx, node_idx);

% get W1, invQ1 and W2
idx2 = 0;
w1k = cell(n_blocks, 1);
invQ1k = cell(n_blocks, 1);
for j = 1 : n_blocks
    idx1 = idx2 + 1;
    idx2 = idx2 + block_size(j);
    w1k{j} = sorted_Q(idx1 : idx2, idx1 : idx2);
    invQ1k{j} = inv(speye(block_size(j)) - c * w1k{j});
end
invQ1 = blkdiag(invQ1k{:});
w2 = sorted_Q - blkdiag(w1k{:});

% do low-rank svds for w2
[u, sigma, v] = svds(w2, r);
v = v';
lambda = inv(inv(sparse(sigma)) - c * v * invQ1 * u);
period_partition_svd = toc(time_partition_svd);

%
% 3. online query 
%
time_online_query = tic;
fprintf(' >> 3. Online Query ... \n');
res = cell(qu_num, 1);
[~, node_idx_inv] = sort(node_idx);
for qi = 1 : qu_num
%     datetime();
    fprintf(' query #%d\n', qi);
    
    s = sparse(qu{qi}, 1, 1 / numel(qu{qi}), n, 1); % select nodes in qu as important nodes
    s = s(node_idx);
    t1 = invQ1 * s;
    t2 = v * t1;
    t3 = lambda * t2;
    t4 = u * t3;
    t5 = invQ1 * t4;
    p = (1 - c) * (t1 + c * t5);
    res{qi} = p(node_idx_inv);
end
period_online_query = toc(time_online_query);

period_total = toc(time_total);

fprintf('\n========== flatApxPPR_MutualExclusion ============\n');

fprintf(' >> Dataset                                 :  %s \n',   ds);
fprintf(' >      # of nodes (certain)                :  %d \n',   n_certain);
fprintf(' >      # of edges (certain)                :  %d \n',   m_certain);
% fprintf(' >      # of possible worlds                :  %d \n',   npw);
fprintf(' >      # of queries                        :  %d \n\n', qu_num);

fprintf(' >> Parameters\n');
fprintf(' >      decay factor (c)                    :  %d \n',   c);
fprintf(' >      low-rank SVD (r)                    :  %d \n',   r);
fprintf(' >      # of partitions                     :  %d \n\n', n_blocks);

fprintf(' >> Total Time                              :  %fs \n', period_total);
fprintf(' >    1) Preprocessing                      :  %f s \n', period_flatten + period_partition_svd);
fprintf(' >         a1) Flatten uncertain graph      :  %f s \n', period_flatten);
fprintf(' >         a2) Partitioning + SVD           :  %f s \n', period_partition_svd);
fprintf(' >    2) Online Query (# of queries = %d)   :  %f s \n\n', qu_num, period_online_query);

fprintf(' >> Average Amortised Total Time per Query              :  %fs \n', period_total / qu_num);
% fprintf(' >> Average Amortised Total Time Time per Query per PW  :  %fs \n', period_total / qu_num / npw);
fprintf(' >> Average Amortised Preprocessing Time per Query      :  %fs \n', (period_flatten + period_partition_svd) / qu_num);


whos
