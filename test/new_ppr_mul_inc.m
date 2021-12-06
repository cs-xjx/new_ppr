function [res] = new_ppr_mul_inc(a, c, kmax, src, tar, qu, ds)
% a = sparse([
% 0, 0, 0, 0, 1, 0, 0, 0, 0;
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
% tar = {[4 7 2 4], [3 5 6 9], [6 9 4], [1 5], [9 2]};
% qu = {[1 2 6]};

time_total = tic;
fprintf('\n >> Start new_ppr_mul_inc ... \n\n');

n_certain = size(a,1); % n_certain: # of nodes in graph
m_certain = nnz(a);

qu_num = numel(qu);
% nsrc = numel(src);

n = max([n_certain, tar{:}]);

% if n > n_certain
    a(n, n) = 0;
% else
%     n = n_certain;
% end

%
% generate possible world (mutual exclusive semantics)
%
time_gen_pwd = tic;
fprintf(' >> Generate Possible World ... \n');
[src_index, tar_nodes, inheritFrom, level, diff_p, isleaf] = tars_inheritfrom_mul(src, tar);
npw = size(inheritFrom, 1);
usrc = unique(src);
nsrc = numel(usrc);
ntar = numel(tar_nodes);
period_gen_pwd = toc(time_gen_pwd);
fprintf('  # of possible worlds            :  %d  \n', npw);
fprintf('  Time (generate possible worlds) :  %fs \n', period_gen_pwd);
%
%

%
% 1. Preprocessing
%
time_preprocessing = tic;
fprintf(' >> 1. Preprocessing ... \n');

% get transition matrix Q and outdegree of src
suma = sum(a, 2);
d = 1 ./ suma;
d(~isfinite(d)) = 0;
Q = a' * spdiags(d, 0, n, n); % transition matrix of a
% D = spdiags(suma(src), 0, nsrc, nsrc);    % outdegree of source nodes
% src_D = suma(src);
% src_D(src_D == 0) = 1;
% D = spdiags(src_D, 0, nsrc, nsrc);

% dispose R and E, add a column of zero vectors to the left-hand side of the matrix
R = inv(eye(n) - c * Q);
% R = [zeros(n, 1) R];
E = speye(n);
% E2 = sparse([1 : n], [2 : n + 1], 1, n, n + 1);
% E2 = [zeros(n, 1) E];

% memorize E(:, tar)*H and E(:, src)*H
sum_H = sparse([], [], 0, n, nsrc, nsrc * nsrc);
sum_H2 = sparse([], [], 0, n, nsrc, nsrc * nsrc);

% save H and pw that we're going to use later on
H_cell = cell(nsrc, 1);
pw_cell = cell(nsrc, 1);

% the treatment of the first possible world (root node)
H_cell{1} = speye(nsrc);
pw_cell{1} = zeros(ntar, 1);

% other possible worlds
for i = 2 : npw
    parent = inheritFrom(i);
    pw_i = pw_cell{level(parent)};
    pw_i(diff_p(i)) = 1;
    insert_i = src_index(diff_p(i));
    idx = find(pw_i);
    
    tar_i = tar_nodes(idx);
    src_id = src_index(idx);
    insert_src = usrc(insert_i);
    M = sparse(tar_i, src_id, 1, n, nsrc);
    L = sparse(src_id, src_id, 1, nsrc, nsrc);
    
    % unique change column
    u = - c * R(usrc, tar_nodes(diff_p(i))) + R(usrc, insert_src);
    if L(insert_i, insert_i) == 1
        u = u + sparse(insert_i, 1, suma(insert_src) - 1, nsrc, 1);
    end
    
    % increment calculation H
    H1 = H_cell{level(parent)};
    y = H1 * u;
    H = H1 - (1 / (1 + y(insert_i))) * y * H1(insert_i, :);
    
    if isleaf(i) == 0
        H_cell{level(i)} = H;
        pw_cell{level(i)} = pw_i;
    end
    
    sum_H = sum_H + sparse(E(:, usrc) * L * H);
    sum_H2 = sum_H2 + sparse(M * H);
    
%     sum_h = sum_h + H * p0(src);
%     sum_h2 = sum_h2 + E(:, tar_i) * H * p0(src);
end
sum_x0 = c * sum_H2 - sum_H;

sum_xk = iter_compute(c, Q, sum_x0, kmax);
FH = (sum_xk + sum_H) / npw;

period_preprocessing = toc(time_preprocessing);

%
% 2. online query 
%
time_online_query = tic;
fprintf(' >> 2. Online Query ... \n');
res = cell(qu_num, 1);
for qi = 1 : qu_num
    s = sparse(qu{qi}, 1, 1 / numel(qu{qi}), n, 1);
    p0 = R * s;
%     p0 = sum(R(:, qu{qi}), 2) / numel(qu{qi});
%     res{qi} = temp_p * p0;
    res{qi} = (1 - c) * (p0 + FH * p0(usrc));
end

period_online_query = toc(time_online_query);

period_total = toc(time_total);


fprintf('\n========== new_ppr_mul_inc ============\n');

fprintf(' >> Dataset                                 :  %s \n',   ds);
fprintf(' >      # of nodes (certain)                :  %d \n',   n_certain);
fprintf(' >      # of edges (certain)                :  %d \n',   m_certain);
fprintf(' >      # of possible worlds                :  %d \n',   npw);
fprintf(' >      # of queries                        :  %d \n\n', qu_num);

fprintf(' >> Parameters\n');
fprintf(' >      decay factor (c)                    :  %d \n',   c);
fprintf(' >      # of iterations                     :  %d \n', kmax);
% fprintf(' >      # of partitions                     :  %d \n', n_blocks);
% fprintf(' >      low-rank SVD (r)                    :  %d \n',   r);
fprintf('\n');

fprintf(' >> Total Time                              :  %fs \n', period_total);
fprintf(' >    1) Preprocessing                      :  %f s \n', period_preprocessing);
fprintf(' >    2) Online Query (# of queries = %d)   :  %f s \n\n', qu_num, period_online_query);

fprintf(' >> Average Amortised Total Time per Query              :  %fs \n', period_total / qu_num);
fprintf(' >> Average Amortised Total Time Time per Query per PW  :  %fs \n', period_total / qu_num / npw);
fprintf(' >> Average Amortised Preprocessing Time per Query      :  %f s \n', period_preprocessing / qu_num);





