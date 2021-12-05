function [res] = UPPR_MutualExclusion(a, c, src, tar, qu, n_blocks, ds) 
time_total = tic;
fprintf('\n >> Start UPPR_MutualExclusion ... \n');

n_certain = size(a, 1); % n_certain: # of nodes in graph
m_certain = nnz(a);

qu_num = numel(qu);
nsrc = numel(src); % nsrc: # of uncertain edge

n = max([n_certain, tar{:}]);

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

tar = cellfun(@(x) [0 x], tar, 'UniformOutput', false);

%
% 1. compute pavg of all pws
%
time_pavg = tic;
fprintf(' >> 1. Compute Pavg of All Possible Worlds ... \n');
% targetSize: # of target nodes of each uncertain edge, store minus without epsilon
targetSize = zeros(1, nsrc);
for i = 1 : nsrc
    tempTarget = tar{i};
    if any(tempTarget == 0)
        targetSize(i) = size(tempTarget, 2);
    else
        targetSize(i) = -size(tempTarget, 2);
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
% compute Pavg
unique_src = unique(src);
unique_src_n = numel(unique_src);
% for i = 1 : n
for src_i = 1 : unique_src_n
    i = unique_src(src_i);
%     c = sum(A(i, :) ~= 0);
    cer_n = floor(sumA(i));
    u = sum(src == i);
%     if u == 0
%         continue;
%     end
    % find n and p of binomial distribution
    ueIndexes = find(src == i);
    ueSize(ueIndexes) = targetSize(ueIndexes);
%     u1 = sum(ueSize ~= 0);
%     ue1Size = ueSize(ueSize ~= 0);
    u1 = sum(ueSize > 0); % u1: # of uncertain edges containing epsilon of a source node
    ue1Size = ueSize(ueSize > 0); % ue1Size: size of each target containing epsilon of a source node
    uniqueSize = unique(ue1Size);
%     unique_n = size(uniqueSize, 2);
%     if unique_n == 1
%         B_n = u1;
%     else
%         B_n = hist(ue1Size, uniqueSize);
%     end
    sorted_ue1Size = sort(ue1Size);
    if numel(sorted_ue1Size) == 0
        sorted_ue1Size = -1;
    end
    B_n = diff(find(diff([-1 sorted_ue1Size sorted_ue1Size(end) + 1])));
    unique_n = numel(B_n);
    
    % X_B: store n and p, one set for each row
    X_B = zeros(size(B_n, 2), 2);
    X_B(:, 1) = B_n';
    X_B(:, 2) = 1 ./ uniqueSize';
    
    u0 = u - u1;
    if cer_n > 0 || u0 > 0
        % compute sum of binomial distribution
        if u1 > 0
            p_h = SumOfBino(X_B, u1);
        else
            p_h = 1;
        end
    
        % same part of certain edges and uncertain edges without epsilon in probability calculation
        Pji = 0;
        for h = 0 : u1
            Pji = Pji + (1 / (cer_n + u - h)) * p_h(h + 1);
        end
        % the probability of certain edges
        certainNode = find(a(i, :) ~= 0);
        P(certainNode, i) = Pji - Tcertain(certainNode, i);
        
        % the probability of uncertain edges without epsilon
        if u0 > 0
            % u0Target: indexes of target of uncertain edges without epsilon
%             u0Target = setdiff(ueIndexes, find(ueSize ~= 0));
            u0Target = find(ueSize < 0);
            for j = 1 : u0
                ue0 = tar{u0Target(j)};
                k = size(ue0, 2);
                ue0 = ue0(ue0 ~= 0);
                P(ue0, i) = P(ue0, i) + Pji / k;  %%%%%
            end
        end
    end
    % the probability of uncertain edges containing epsilon
    for j = 1 : unique_n
        % tempX_B: store n and p of different uncertain edges containing epsilon
        tempX_B = X_B;
        if B_n(j) == 1
            tempX_B(j, :) = [];
        else
            tempX_B(j, 1) = B_n(j) - 1;
        end
        if u1 - 1 == 0
            p_h = 1;
        else
            p_h = SumOfBino(tempX_B, u1 - 1);
        end
        Pji = 0;
        for h = 0 : u1 - 1
            Pji = Pji + (1 / (cer_n + u - h)) * p_h(h + 1);
        end
        k = uniqueSize(j);
        u1Target = find(ueSize == k);
        for l = 1 : B_n(j)
            ue1 = tar{u1Target(l)};
            ue1 = ue1(ue1 ~= 0);
            P(ue1, i) = P(ue1, i) + Pji / k; %%%%%
        end
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
    
    uppr_mut = (1-c)*t1 + c*(1-c)*t3 + (c*(1-c))*t5 + c*c*(1-c)*t7 + (c*c*(1-c))*t9 + (c*c*(1-c))*t11;
    
    res{qi} = uppr_mut(node_idx_inv);
%     [x, y] = sort(uppr_mut, 'descend');
%     x(1 : 10)
%     y(1 : 10)
%     full(uppr_mut)
end
period_online_query = toc(time_online_query);

period_total = toc(time_total);

fprintf('\n========== UPPR_MutualExclusion ============\n');

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
