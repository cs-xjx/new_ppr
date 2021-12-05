function [res] = exhApxPPR_reference_MutualExclusion(a, c, n_blocks, src, tar, qu, r, ds) 
time_total = tic;
fprintf('\n >> Start exhApxPPR2_MutualExclusion ... \n');

% fpath = 'D:\Polyspace\metis-5.1.0\metismex-dgleich\dataset\';
% fname = 'ego-Facebook';
% fn = [fpath, fname, '.mat'];
% 
% load(fn);
% A = Problem.A;
% src = [1912 3059];
% tar = {[3356 2221 1154], [1534 306 4071]};
% qu = {[260 553], [2658 145], [1163 3525], [2930 3248], [1457 1100]};

% Problem = load('digraph_random_n100_e200.mat');
% A = Problem.A;
% source = [76, 22, 18, 49];
% target = {[-1, 51, 6, 43], [-1, 23, 15, 19], [-1, 31, 45, 33], [-1, 2, 43, 46]};
% s = sparse([90 94 22 28 41 18], 1, 1/6, 100, 1);

% A = sparse([
%     0, 0, 0, 0, 1, 0, 0, 0, 0;
% 0, 0, 1, 1, 0, 0, 0, 0, 0;
% 0, 0, 0, 0, 0, 1, 0, 0, 0;
% 0, 0, 0, 0, 1, 0, 0, 0, 1;
% 0, 1, 0, 0, 0, 0, 1, 0, 0;
% 0, 1, 0, 0, 0, 0, 0, 0, 0;
% 0, 0, 1, 0, 0, 0, 0, 0, 0;
% 0, 0, 0, 1, 0, 0, 0, 0, 0;
% 1, 0, 0, 0, 0, 0, 0, 0, 0
%     ]);
% src = [1 2 1 2 5];
% tar = {[2 7 14 17], [119 8 6 17], [4 9 6], [90 5], [9 117]};
% qu = {[1 2 6]};

% A = sparse([
% 0 1 0 0 1;
% 0 0 1 1 0;
% 0 0 0 0 1;
% 0 0 1 0 0;
% 0 0 0 0 0
%     ]);
% src = [1,4,5];
% tar = {[2,3,5],[3,4],[2 5]};
% qu = {[3 5]};
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % common parameter 
% c = 0.8;            % decay factor
% kmax = 20;         % # of maximum iterations
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

n_certain = size(a, 1);
m_certain = nnz(a);

qu_num = numel(qu);
nsrc = numel(src); 

% numUE: # of uncertain edge
% tar = cellfun(@(x) [0 x], tar, 'UniformOutput', false); % if target nodes include epsilon, the fun runs, and zero represrnts epsilon
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% begin = tic;

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
tars = cat(nsrc+1, tars{:});
tars = reshape(tars,[],nsrc);  % tars: possible world
npw = size(tars, 1);
period_gen_pwd = toc(time_gen_pwd);
fprintf('  # of possible worlds            :  %d  \n', npw);
fprintf('  Time (generate possible worlds) :  %fs \n', period_gen_pwd);

fprintf('\n== BEGIN online query ==\n');

% initilising
res = cell(qu_num, 1);
for qi = 1:qu_num
    res{qi} = sparse(n, 1);
end

% enumerate all possible worlds
period_partition = 0;
period_svd_lambda = 0;
period_online_query = 0;

for i = 1:npw
    
    time_partition = tic;
    
    % find uncertain part of pw
    tar_i = tars(i,:);
    idx = find(tar_i);
    tar_i = tar_i(idx);
    src_i = src(idx);

    % get pq (transition matrix of pw)
    pa = a + sparse(src_i, tar_i, 1, n, n);
    d = 1./sum(pa,2);
    d(~isfinite(d)) = 0;
    pq = pa' * spdiags(d,0,n,n);
    clear d pa;

    % partition pq
    [map, ~] = metismex('PartGraphRecursive', pq, n_blocks); % small n_blocks
    %[map, ~] = metismex('PartGraphKway', T, k);             % large n_blocks
    [sorted_map, node_idx] = sort(map);
    block_size = diff(find(diff([-1 sorted_map n])));

    % orderDiff = diff(nodeOrder);
    % ltzeroIndex = find(orderDiff < 0);
    % kNum = diff(ltzeroIndex);
    % kNum = [ltzeroIndex(1) kNum n - ltzeroIndex(end)];

    sorted_pq = pq(node_idx, node_idx);

    % tic
%     t_invQ1k = cell(1, n_blocks);
%     blk_reorderT = mat2cell(sorted_pq, block_size, block_size);
%     for j = 1 : n_blocks
%         blk_size = block_size(j);
%         t_invQ1k{j} = inv(speye(blk_size) - c * blk_reorderT{j,j});
%         blk_reorderT{j,j} = sparse(blk_size, blk_size);
%     end
%     invQ1 = blkdiag(t_invQ1k{:});
%     w2 = cell2mat(blk_reorderT);
    % toc

    % get W1, invQ1 and W2

    idx2 = 0;
    w1 = cell(1, n_blocks);
    invQ1k = cell(1, n_blocks);
    for j = 1 : n_blocks
        idx1 = idx2 + 1;
        idx2 = idx2 + block_size(j);
        w1{j} = sorted_pq(idx1:idx2, idx1:idx2);
        invQ1k{j} = inv(speye(block_size(j)) - c * w1{j});
    end    
    invQ1 = blkdiag(invQ1k{:});
    w2 = sorted_pq - blkdiag(w1{:});
    period_partition = period_partition + toc(time_partition);

    % low-rank svd for w2
    time_svd_lambda = tic;
    [u, sigma, v] = svds(w2, r);
    v = v';
    lambda = inv(inv(sparse(sigma)) - c*v*invQ1*u);
    period_svd_lambda = period_svd_lambda + toc(time_svd_lambda);

    %[temp_period_partAndInv, temp_period_svds, temp_period_A, invQ1, v, lambda, u, node_idx] = ApxPPR_B_LIN(pq, c, n_blocks, r);

    % online query
    time_online_query = tic;   
    [~, node_idx_inv] = sort(node_idx);   
    for qi = 1 : qu_num
        %fprintf(' query = #%d,  pwd = #%d\n', qi, i);
        s = sparse(qu{qi}, 1, 1 / numel(qu{qi}), n, 1); % select nodes in qu as important nodes
        s = s(node_idx);

        qs = invQ1*s;       
        vqs = v * qs;
        avqs = lambda * vqs;
        uavqs = u * avqs;
        quavqs = invQ1 * uavqs;
        psum = (1 - c) * (qs + c * quavqs);

%        psum = (1-c)*(qs + c*invQ1*u*lambda*v*qs);
        psum = psum(node_idx_inv);

        %%%  query
        res{qi} = res{qi} + psum;
    end
    period_online_query = period_online_query + toc(time_online_query);
    
% % % % % % % % % % % % % % % % %         
%        pSum = pSum + temp_pSum;
%        period_partAndInv = period_partAndInv + temp_period_partAndInv;
%        period_svds = period_svds + temp_period_svds;
%        period_a = period_a + temp_period_A;
end
%     p_eAP_Mut = pSum / npw;
%     res{qi} = pSum / npw;

for qi = 1 : qu_num
    res{qi} = res{qi} / npw;
end
fprintf('== END online query ==\n');
period_total = toc(time_total);

fprintf('\n========== exhApxPPR2_MutualExclusion ============\n');

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
fprintf(' >    1) Preprocessing                      :  %fs \n', period_partition+period_svd_lambda);
fprintf(' >         a1) Partitioning & Inverse       :  %fs \n', period_partition);
fprintf(' >         a2) SVD & Lambda = S-cVQU        :  %fs \n', period_svd_lambda);
fprintf(' >    2) Online Query (# of queries = %d)   :  %fs \n\n', qu_num, period_online_query);

fprintf(' >> Average Amortised Total Time per Query              :  %fs \n', period_total / qu_num);
fprintf(' >> Average Amortised Total Time Time per Query per PW  :  %fs \n', period_total / qu_num / npw);
fprintf(' >> Average Amortised Preprocessing Time per Query      :  %fs \n', (period_partition+period_svd_lambda) / qu_num);

