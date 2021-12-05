function [res] = collApx2PPR_MutualExclusion(a, c, n_blocks, src, tar, qu, r, ds) 
time_total = tic;
fprintf('\n >> Start collApx2PPR_MutualExclusion ... \n');

n_certain = size(a, 1);
m_certain = nnz(a);

qu_num = numel(qu);
nsrc = numel(src); % numUE: # of uncertain edges

n = max([n_certain, tar{:}]);
a(n, n) = 0;

%
% generate possible world (mutual exclusive semantics)
%
time_gen_pwd = tic;
fprintf(' >> Generate Possible World ... \n');
tar0 = cellfun(@(x) [0 x], tar, 'UniformOutput', false);
tars = cell(1,nsrc);
[tars{end:-1:1}] = ndgrid(tar0{end:-1:1});
tars = cat(nsrc + 1, tars{:});
tars = reshape(tars,[],nsrc);  % tars: possible world
npw = size(tars, 1);
period_gen_pwd = toc(time_gen_pwd);
fprintf('  # of possible worlds            :  %d  \n', npw);
fprintf('  Time (generate possible worlds) :  %fs \n', period_gen_pwd);

%
% 1. get w1_ave and w2_ave 
%
time_average_w1_w2 = tic;
fprintf(' >> 1. Get w1_ave and w2_ave ... \n');
w1_ave = sparse(n, n);
w2_ave = sparse(n, n);
for i = 1:npw
    % find the uncertain part of the pw
    tar_i = tars(i,:);
    idx = find(tar_i);
    tar_i = tar_i(idx);
    src_i = src(idx);

    % pq: transition matrix
    pa = a + sparse(src_i, tar_i, 1, n, n);
    d = 1./sum(pa,2);
    d(~isfinite(d)) = 0;
    pq = pa' * spdiags(d,0,n,n);
    clear d pa;

    % get w1 and w2
    % partition pq
%     [map, ~] = metismex('PartGraphRecursive', pq, n_blocks); % small n_blocks
    [map, ~] = metismex('PartGraphKway', pq, n_blocks);             % large n_blocks
    [sorted_map, node_idx] = sort(map);
    block_size = diff(find(diff([-1 sorted_map n])));
    sorted_pq = pq(node_idx, node_idx);    
    idx2 = 0;
    w1_blks = cell(1, n_blocks);
    for j = 1:n_blocks
        idx1 = idx2 + 1;
        idx2 = idx2 + block_size(j);
        w1_blks{j} = sorted_pq(idx1:idx2, idx1:idx2);
    end
    w1 = blkdiag(w1_blks{:});
    w2 = sorted_pq - w1;
    [~, node_idx_inv] = sort(node_idx);
    w1 = w1(node_idx_inv, node_idx_inv);
    w2 = w2(node_idx_inv, node_idx_inv);
    
    w1_ave = w1_ave + w1;
    w2_ave = w2_ave + w2;
end
w1_ave = w1_ave / npw;
w2_ave = w2_ave / npw;
% do low-rank svds for w2
[u, sigma, v] = svds(w2_ave, r);
v = v';
invQ1 = inv(speye(n) - c * w1_ave);
lambda = inv(inv(sparse(sigma)) - c * v * invQ1 * u);
period_average_w1_w2 = toc(time_average_w1_w2);

%
% 2. online query 
%
time_online_query = tic;
fprintf(' >> 2. Online Query ... \n');
res = cell(qu_num, 1);
for qi = 1 : qu_num
    fprintf(' query #%d\n', qi);
    s = sparse(qu{qi}, 1, 1 / numel(qu{qi}), n, 1); % select nodes in qu as important nodes
    t1 = invQ1 * s;
    t2 = v * t1;
    t3 = lambda * t2;
    t4 = u * t3;
    t5 = invQ1 * t4;
    res{qi} = (1 - c) * (t1 + c * t5);
end
period_online_query = toc(time_online_query);
period_total = toc(time_total);

fprintf('\n========== collApx2PPR_MutualExclusion ============\n');

fprintf(' >> Dataset                                 :  %s \n',   ds);
fprintf(' >      # of nodes (certain)                :  %d \n',   n_certain);
fprintf(' >      # of edges (certain)                :  %d \n',   m_certain);
fprintf(' >      # of possible worlds                :  %d \n',   npw);
fprintf(' >      # of queries                        :  %d \n\n', qu_num);

fprintf(' >> Parameters\n');
fprintf(' >      decay factor (c)                    :  %d \n',   c);
fprintf(' >      low-rank SVD (r)                    :  %d \n',   r);
fprintf(' >      # of partitions                     :  %d \n\n', n_blocks);

fprintf(' >> Total Time                                        :  %f s \n', period_total);
fprintf(' >    1) Preprocessing (Ave w1, w2 + Partition + SVD) :  %f s \n', period_average_w1_w2);
fprintf(' >    2) Online Query (# of queries = %d)             :  %f s \n\n', qu_num, period_online_query);

fprintf(' >> Average Amortised Total Time per Query            :  %fs \n', period_total / qu_num);
fprintf(' >> Average Amortised Total Time per Query per PW     :  %fs \n', period_total / qu_num / npw);
fprintf(' >> Average Amortised Preprocessing Time per Query    :  %fs \n', period_average_w1_w2 / qu_num);


whos
