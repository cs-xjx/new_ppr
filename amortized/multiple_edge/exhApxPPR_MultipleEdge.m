function [res] = exhApxPPR_MultipleEdge(a, c, n_blocks, src, tar, qu, r, ds) 
time_total = tic;
fprintf('\n >> Start exhApxPPR_MultipleEdge ... \n');

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
% tar0 = cellfun(@(x) [0 x], tar, 'UniformOutput', false);
tars = cartprod_MultipleEdge(tar, nsrc);
[npw, ~] = size(tars);
period_gen_pwd = toc(time_gen_pwd);
fprintf('  # of possible worlds            :  %d  \n', npw);
fprintf('  Time (generate possible worlds) :  %fs \n', period_gen_pwd);
%
%

%
% 1. partition and online query
%
fprintf(' >> 1. Partition and Online Query ... \n');
res = cell(qu_num, 1);
for qi = 1 : qu_num
    res{qi} = sparse(n, 1);
end

% enumerate all possible worlds by loop
period_partition = 0;
period_svd_lambda = 0;
period_online_query = 0;

for i = 1 : npw
    
    time_partition = tic;
    
    pa = a;
    for j = 1 : nsrc
        pa = pa + sparse(src(j), tars{i, j}, 1, n, n);
    end
    % pq: transition matrix
    d = 1./sum(pa,2);
    d(~isfinite(d)) = 0;
    pq = pa' * spdiags(d,0,n,n);
    clear d pa;
    
    % partition pq
    if n_blocks > n
        error('the value of k should not greater than n');
    end
%     [map, ~] = metismex('PartGraphRecursive', pq, n_blocks); % small n_blocks
    [map, ~] = metismex('PartGraphKway', pq, n_blocks);     % large n_blocks
    [sorted_map, node_idx] = sort(map);
    block_size = diff(find(diff([-1 sorted_map n])));
    
    sorted_pq = pq(node_idx, node_idx);
    
    % get W1, invQ1 and W2
    idx2 = 0;
    w1k = cell(n_blocks, 1);
    invQ1k = cell(n_blocks, 1);
    for j = 1 : n_blocks
        idx1 = idx2 + 1;
        idx2 = idx2 + block_size(j);
        w1k{j} = sorted_pq(idx1 : idx2, idx1 : idx2);
        invQ1k{j} = inv(speye(block_size(j)) - c * w1k{j});
    end
    invQ1 = blkdiag(invQ1k{:});
    w2 = sorted_pq - blkdiag(w1k{:});
    clear invQ1k w1k pq;
    % invQ1
    period_partition = period_partition + toc(time_partition);
    
    time_svd_lambda = tic;
    % do low-rank svds for w2
    [u, sigma, v] = svds(w2, r);
    clear w2;
    v = v';
    
    lambda = inv(inv(sparse(sigma)) - c * v * invQ1 * u);
    period_svd_lambda = period_svd_lambda + toc(time_svd_lambda);
%     [temp_period_partAndInv, temp_period_svds, temp_period_A, invQ1, v, A, u, nodeOrder] = ApxPPR_B_LIN(pq, c, n_blocks, r);    
    % online query
    time_online_query = tic;
    [~, node_idx_inv] = sort(node_idx);
    for qi = 1 : qu_num
        s = sparse(qu{qi}, 1, 1 / numel(qu{qi}), n, 1); % select nodes in qu as important nodes
        s = s(node_idx);
        
        t1 = invQ1 * s;
        t2 = v * t1;
        t3 = lambda * t2;
        t4 = u * t3;
        t5 = invQ1 * t4;
        psum = (1 - c) * (t1 + c * t5);
        psum = psum(node_idx_inv);
        res{qi} = res{qi} + psum;
%         res{qi} = res{qi} + B_LIN_ppr(c, invQ1, v, A, u, nodeOrder, nodeOrder2, qu{qi});
    end
    period_online_query = period_online_query + toc(time_online_query);
    clear invQ1 u v;
end

time_online_query = tic;
for qi = 1 : qu_num
    res{qi} = res{qi} / npw;
end
period_online_query = period_online_query + toc(time_online_query);

period_total = toc(time_total);

fprintf('\n========== exhApxPPR_MultipleEdge ============\n');

fprintf(' >> Dataset                                 :  %s \n',   ds);
fprintf(' >      # of nodes (certain)                :  %d \n',   n_certain);
fprintf(' >      # of edges (certain)                :  %d \n',   m_certain);
fprintf(' >      # of possible worlds                :  %d \n',   npw);
fprintf(' >      # of queries                        :  %d \n\n', qu_num);

fprintf(' >> Parameters\n');
fprintf(' >      decay factor (c)                    :  %d \n',   c);
fprintf(' >      low-rank SVD (r)                    :  %d \n',   r);
fprintf(' >      # of partitions                     :  %d \n\n', n_blocks);

fprintf(' >> Total Time                              :  %fs \n', period_total);
fprintf(' >    1) Preprocessing                      :  %fs \n', period_partition + period_svd_lambda);
fprintf(' >         a1) Partitioning & Inverse       :  %fs \n', period_partition);
fprintf(' >         a2) SVD & Lambda = S-cVQU        :  %fs \n', period_svd_lambda);
fprintf(' >    2) Online Query (# of queries = %d)   :  %fs \n\n', qu_num, period_online_query);

fprintf(' >> Average Amortised Total Time per Query              :  %fs \n', period_total / qu_num);
fprintf(' >> Average Amortised Total Time Time per Query per PW  :  %fs \n', period_total / qu_num / npw);
fprintf(' >> Average Amortised Preprocessing Time per Query      :  %fs \n', (period_partition + period_svd_lambda) / qu_num);


whos


