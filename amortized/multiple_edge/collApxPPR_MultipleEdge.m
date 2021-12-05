function [res] = collApxPPR_MultipleEdge(a, c, n_blocks, src, tar, qu, r, ds) 
time_total = tic;
fprintf('\n >> Start collApxPPR_MultipleEdge ... \n');

n_certain = size(a, 1);
m_certain = nnz(a);

qu_num = numel(qu);
nsrc = numel(src); % numUE: # of uncertain edge

n = max([n_certain, tar{:}]);
a(n, n) = 0;

%
% generate possible world (multiple edge semantics)
%
time_gen_pwd = tic;
fprintf(' >> Generate Possible World ... \n');
tars = cartprod_MultipleEdge(tar, nsrc);
[npw, ~] = size(tars);
period_gen_pwd = toc(time_gen_pwd);
fprintf('  # of possible worlds            :  %d  \n', npw);
fprintf('  Time (generate possible worlds) :  %fs \n', period_gen_pwd);
%
%

%
% 1. collapse all pws
%
time_collapse = tic;
fprintf(' >> 1. Collapse All Possible Worlds ... \n');
qavg = sparse(n, n);
% enumerate all possible worlds by loop
for i = 1 : npw
    pa = a;
    for j = 1 : nsrc
        pa = pa + sparse(src(j), tars{i, j}, 1, n, n);
    end
    % pq: transition matrix
    d = 1./sum(pa,2);
    d(~isfinite(d)) = 0;
    pq = pa' * spdiags(d,0,n,n);
    clear d pa;
    % get the avg of transition matrixes of all pws
    qavg = qavg + pq / npw;
end
period_collapse = toc(time_collapse);

%
% 2. partition collapsed PW (qavg)
%
time_partition_svd = tic;
fprintf(' >> 2. Partition Collapsed PW ... \n');
% [map, ~] = metismex('PartGraphRecursive', qavg, n_blocks); % small n_blocks
[map, ~] = metismex('PartGraphKway', qavg, n_blocks);     % large n_blocks
[sorted_map, node_idx] = sort(map);
block_size = diff(find(diff([-1 sorted_map n])));
sorted_pq_avg = qavg(node_idx, node_idx);

% get W1, invQ1 and W2
idx2 = 0;
w1k = cell(n_blocks, 1);
invQ1k = cell(n_blocks, 1);
for j = 1 : n_blocks
    idx1 = idx2 + 1;
    idx2 = idx2 + block_size(j);
    w1k{j} = sorted_pq_avg(idx1 : idx2, idx1 : idx2);
    invQ1k{j} = inv(speye(block_size(j)) - c * w1k{j});
end
invQ1 = blkdiag(invQ1k{:});
w2 = sorted_pq_avg - blkdiag(w1k{:});

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
    fprintf('     query #%d\n', qi);
    
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

fprintf('\n========== collApxPPR_MultipleEdge ============\n');

fprintf(' >> Dataset                                 :  %s \n',   ds);
fprintf(' >      # of nodes (certain)                :  %d \n',   n_certain);
fprintf(' >      # of edges (certain)                :  %d \n',   m_certain);
fprintf(' >      # of possible worlds                :  %d \n',   npw);
fprintf(' >      # of queries                        :  %d \n\n', qu_num);

fprintf(' >> Parameters\n');
fprintf(' >      decay factor (c)                    :  %d \n',   c);
fprintf(' >      low-rank SVD (r)                    :  %d \n',   r);
fprintf(' >      # of partitions                     :  %d \n\n', n_blocks);

fprintf(' >> Total Time                              :  %f s \n', period_total);
fprintf(' >    1) Preprocessing                      :  %f s \n', period_collapse + period_partition_svd);
fprintf(' >         a1) Collapse All Possible Worlds :  %f s \n', period_collapse);
fprintf(' >         a2) Partitioning + SVD           :  %f s \n', period_partition_svd);
fprintf(' >    2) Online Query (# of queries = %d)   :  %f s \n\n', qu_num, period_online_query);

fprintf(' >> Average Amortised Total Time per Query              :  %f s \n', period_total / qu_num);
fprintf(' >> Average Amortised Total Time per Query per PW       :  %f s \n', period_total / qu_num / npw);
fprintf(' >> Average Amortised Preprocessing Time per Query      :  %f s \n', (period_collapse + period_partition_svd) / qu_num);

whos

