function [res] = UPPR_MultipleEdge(a, c, src, tar, qu, n_blocks, ds) 
time_total = tic;
fprintf('\n >> Start UPPR_MultipleEdge ... \n');

n_certain = size(a, 1); % n_certain: # of nodes in graph
m_certain = nnz(a);

qu_num = numel(qu);
nsrc = numel(src); % nsrc: # of uncertain edge

n = max([n_certain, tar{:}]);

a(n, n) = 0;

% inserts the newly added node into the graph
% new_nodes = cellfun(@(x) x(x > n), tar, 'UniformOutput', false);
% new_nodes = cat(2, new_nodes{:});
% new_nodes_num = 0;
% while numel(new_nodes) ~= 0
%     node = min(new_nodes);
%     new_nodes_num = new_nodes_num + 1;
%     for i = 1 : numel(tar)
%         tar{i}(tar{i} == node) = n + new_nodes_num;
%     end
%     new_nodes = new_nodes(new_nodes ~= node);
% end
% if new_nodes_num > 0
%     A = blkdiag(A, zeros(new_nodes_num));
%     n = n + new_nodes_num;
% end

tar = cellfun(@(x) [0 x], tar, 'UniformOutput', false);

%
% 1. compute pavg of all pws
%
time_pavg = tic;
fprintf(' >> 1. Compute Pavg of All Possible Worlds ... \n');
% targetSize: # of target nodes of each uncertain edge, not including epsilon nodes
targetSize = zeros(1, nsrc);
for i = 1 : nsrc
    tempTarget = tar{i};
    if any(tempTarget == 0)
        targetSize(i) = size(tempTarget, 2) - 1;
    else
        targetSize(i) = size(tempTarget, 2);
    end
end

% Tcertain: certain parts of the graph
% cer_n: # of certain edges of a source node
% u: # of uncertain edges of a source node
sumA = sum(a, 2);
d = 1 ./ sumA;
d(~isfinite(d)) = 0;
Tcertain = a' * spdiags(d, 0, n, n);
clear d;

% ueSize: size of each target of a source node
ueSize = zeros(1, nsrc);
P = sparse(n, n); % P: Pavg of UPPR
unique_src = unique(src);
unique_src_n = numel(unique_src);
% compute Pavg
for src_i = 1 : unique_src_n
    i = unique_src(src_i);
    cer_n = floor(sumA(i));
%     u = sum(src == i);
%     if u == 0
%         continue;
%     end
    % find n and p of binomial distribution
    ueIndexes = find(src == i);
    ueSize(ueIndexes) = targetSize(ueIndexes);
    
    % X_B: store n and p, one set for each row
    % B_n: # of target nodes of uncertain edges
    u1 = sum(ueSize);
    X_B = zeros(1, 2);
    X_B(1, 2) = 1 / 2;
    
    % the probability of certain edges
    if cer_n > 0
        X_B(1, 1) = u1;
        % compute sum of binomial distribution
        p_h = SumOfBino(X_B, u1);
        
        % probability calculation of choose epsilon
        Pji = 0;
        for h = 0 : u1
            Pji = Pji + (1 / (cer_n + h)) * p_h(h + 1);
        end
        % the probability of certain edges
        certainNode = find(a(i, :) ~= 0);
        P(certainNode, i) = Pji - Tcertain(certainNode, i);
    end

    % the probability of uncertain edges
    X_B(1, 1) = u1 - 1;
    p_h = SumOfBino(X_B, X_B(1, 1));
    Pji = 0;
    for h = 0 : X_B(1, 1)
        Pji = Pji + (1 / (cer_n + 1 + h)) * p_h(h + 1);
    end
    tar_num = numel(ueIndexes);
    for j = 1 : tar_num
        ue = tar{ueIndexes(j)};
        ue = ue(ue ~= 0);
        P(ue, i) = P(ue, i) + Pji / 2;
    end
    ueSize(:) = 0;
end
period_pavg = toc(time_pavg);

%
% 2. partition Tcertain
%
time_partition_svd = tic;
fprintf(' >> 2. Partition Tcertain ... \n');
if n_blocks > n
    error('the value of k should not greater than n');
end
% [map, ~] = metismex('PartGraphRecursive', sparse(Tcertain), n_blocks); % small n_blocks
[map, ~] = metismex('PartGraphKway', sparse(Tcertain), n_blocks); % large n_blocks

[sorted_map, node_idx] = sort(map);
block_size = diff(find(diff([-1 sorted_map n])));
sorted_T = Tcertain(node_idx, node_idx);

% get TBL, invQBL and TX
idx2 = 0;
TBLk = cell(n_blocks, 1);
invQBLk = cell(n_blocks, 1);
for i = 1 : n_blocks
    idx1 = idx2 + 1;
    idx2 = idx2 + block_size(i);
    TBLk{i} = sorted_T(idx1 : idx2, idx1 : idx2);
    invQBLk{i} = inv(speye(block_size(i)) - c * TBLk{i});
end
invQBL = blkdiag(invQBLk{:});
TX = sorted_T - blkdiag(TBLk{:});
period_partition_svd = toc(time_partition_svd);

%
% 3. online query 
%
time_online_query = tic;
fprintf(' >> 3. Online Query ... \n');
P = P(node_idx, node_idx);
res = cell(qu_num, 1);
[~, node_idx_inv] = sort(node_idx);
for qi = 1 : qu_num
    fprintf(' query #%d\n', qi);
%     [uppr_mut, temp_period_partAndInv] = metis_uppr(Tcertain, P, c, n_blocks, qu{qi});
    s = sparse(qu{qi}, 1, 1 / numel(qu{qi}), n, 1); % select nodes in qu as important nodes
    s = s(node_idx);

    t1 = invQBL * s; %1
    t2 = TX * t1; %2
    t3 = invQBL * t2; %2
    t4 = P * t1; %3
    t5 = invQBL * t4; %3
    t6 = TX * t3; %4
    t7 = invQBL * t6; %4
    t8 = P * t3; %5
    t9 = invQBL * t8; %5
    t10 = TX * t5; %6
    t11 = invQBL * t10; %6
    
    uppr_mul = (1-c)*t1 + c*(1-c)*t3 + (c*(1-c))*t5 + c*c*(1-c)*t7 + (c*c*(1-c))*t9 + (c*c*(1-c))*t11;
    
    res{qi} = uppr_mul(node_idx_inv);
end
period_online_query = toc(time_online_query);

period_total = toc(time_total);

fprintf('\n========== UPPR_MultipleEdge ============\n');

fprintf(' >> Dataset                                 :  %s \n',   ds);
fprintf(' >      # of nodes (certain)                :  %d \n',   n_certain);
fprintf(' >      # of edges (certain)                :  %d \n',   m_certain);
% fprintf(' >      # of possible worlds                :  %d \n',   npw);
fprintf(' >      # of queries                        :  %d \n\n', qu_num);

fprintf(' >> Parameters\n');
fprintf(' >      decay factor (c)                    :  %d \n',   c);
% fprintf(' >      low-rank SVD (r)                    :  %d \n',   r);
fprintf(' >      # of partitions                     :  %d \n\n', n_blocks);

fprintf(' >> Total Time                              :  %fs \n', period_total);
fprintf(' >    1) Preprocessing                      :  %f s \n', period_pavg + period_partition_svd);
fprintf(' >         a1) Compute pavg                 :  %f s \n', period_pavg);
fprintf(' >         a2) Partitioning + SVD           :  %f s \n', period_partition_svd);
fprintf(' >    2) Online Query (# of queries = %d)   :  %f s \n\n', qu_num, period_online_query);

fprintf(' >> Average Amortised Total Time per Query              :  %fs \n', period_total / qu_num);
fprintf(' >> Average Amortised Preprocessing Time per Query      :  %fs \n', (period_pavg + period_partition_svd) / qu_num);

whos
