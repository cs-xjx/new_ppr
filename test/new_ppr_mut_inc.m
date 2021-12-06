function [res] = new_ppr_mut_inc(a, c, kmax, src, tar, qu, ds)

time_total = tic;
fprintf('\n >> Start new_ppr_mut_inc ... \n\n');

n_certain = size(a,1); % n_certain: # of nodes in graph
m_certain = nnz(a);

qu_num = numel(qu);
nsrc = numel(src);

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
% [src, sorted_idx] = sort(src);
% tar = tar(sorted_idx);
tar0 = cellfun(@(x) [0 x], tar, 'UniformOutput', false);
% tar0 = tar;
tars = cell(1, nsrc);
[tars{end:-1:1}] = ndgrid(tar0{end:-1:1});
tars = cat(nsrc + 1, tars{:});
tars = reshape(tars, [], nsrc);  % tars: possible world
npw = size(tars, 1);
[inheritFrom, level, diff_p, isleaf] = inheritance_relation(tar0, npw);
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
src_D = suma(src);
src_D(src_D == 0) = 1;
% D = spdiags(src_D, 0, nsrc, nsrc);

% dispose R and E, add a column of zero vectors to the left-hand side of the matrix
R = inv(speye(n) - c * Q);
% R = [zeros(n, 1) R];
E = speye(n);
% E2 = sparse([1 : n], [2 : n + 1], 1, n, n + 1);
E2 = [zeros(n, 1) E];

% memorize E(:, tar)*H and E(:, src)*H
sum_H = sparse([], [], 0, n, nsrc, nsrc * nsrc);
sum_H2 = sparse([], [], 0, n, nsrc, nsrc * nsrc);

% s = sparse(qu{1}, 1, 1 / numel(qu{1}), n, 1);
% p0 = R(:, 2 : n + 1) * s;
% sum_h = zeros(nsrc, 1);
% sum_h2 = sparse(n, 1);

% save H that we're going to use later on
H_cell = cell(nsrc, 1);
% h_cell = cell(nsrc, 1);

% the treatment of the first possible world (root node)
tar_i = tars(1, :);
idx_zero = tar_i == 0;
idx = ~idx_zero;
src_i = src;
src_i(idx_zero) = 0;
% src_D(intersect(find(src_D == 0), idx_zero)) = 1;
% D = spdiags(src_D, 0, nsrc, nsrc);
D = diag(src_D);

D(:, idx) = D(:, idx) - c * R(src, tar_i(idx)) + R(src, src_i(idx));
H = inv(D);
% H = inv(D - c * R(src, tar_i + 1) + R(src, src_i + 1));
% h = H * E(src, :);
H_cell{1} = H;
% h_cell{1} = h;
% clear H;
% sum_h = sum_h + h;
sum_H = sum_H + sparse(E2(:, src_i + 1) * H);
sum_H2 = sum_H2 + sparse(E2(:, tar_i + 1) * H);
% sum_h = sum_h + H * p0(src);
% sum_h2 = sum_h2 + E(:, tar_i) * H * p0(src);

% other possible worlds
for i = 2 : npw
    tar_i = tars(i, :);
    parent = inheritFrom(i);
    insert_i = diff_p(i);
    if parent <= inheritFrom(i - 1)
        idx_zero = tar_i == 0;
        src_i = src;
        src_i(idx_zero) = 0;
    else
        src_i(insert_i) = src(insert_i);
    end
    % unique change column
%     u = sparse(insert_i, 1, suma(src_i(insert_i)) - src_D(insert_i), nsrc, 1) - c * (R(src, tar_i(insert_i) + 1) - R(src, tars(parent, insert_i) + 1));
    u = - c * R(src, tar_i(insert_i));
    u(insert_i) = u(insert_i) + suma(src_i(insert_i)) - src_D(insert_i);
%     u = R(src, tar_i(insert_i) + 1) - R(src, tars(parent, insert_i) + 1);
    if tars(parent, insert_i) == 0
        u = u + R(src, src_i(insert_i));
    else
        u = u + c * R(src, tars(parent, insert_i));
    end
    
    % increment calculation H
    H1 = H_cell{level(parent)};
    y = H1 * u;
    H = H1 - (1 / (1 + y(insert_i))) * y * H1(insert_i, :);
%     H = H1 + c / (1 - c * y(insert_i)) * y * H1(insert_i, :);
    
    if isleaf(i) == 0
        H_cell{level(i)} = H;
    end
    
    sum_H = sum_H + sparse(E2(:, src_i + 1) * H);
    sum_H2 = sum_H2 + sparse(E2(:, tar_i + 1) * H);
    
%     sum_h = sum_h + H * p0(src);
%     sum_h2 = sum_h2 + E(:, tar_i) * H * p0(src);
end
sum_x0 = c * sum_H2 - sum_H;

sum_xk = iter_compute(c, Q, sum_x0, kmax);
clear sum_x0 Q;
FH = (sum_xk + sum_H) / npw;

% R(:, 1) = [];

period_preprocessing = toc(time_preprocessing);

%
% 2. online query 
%
time_online_query = tic;
fprintf(' >> 2. Online Query ... \n');
res = cell(qu_num, 1);
for qi = 1 : qu_num
    s = sparse(qu{qi}, 1, 1 / numel(qu{qi}), n, 1);
%     p0 = R(:, 2 : n + 1) * s;
    p0 = sum(R(:, qu{qi}), 2) / numel(qu{qi});
%     res{qi} = temp_p * p0;
    res{qi} = (1 - c) * (p0 + FH * p0(src));
end

% I_sh = E(:, src) * sum_h;
% sum_x0 = c * sum_h2 - I_sh;
% sum_xk = iter_compute(c, Q, sum_x0, kmax);
% res{1} = p0 + sum_xk / npw + I_sh / npw;
% res{1} = (1 - c) * res{1};

period_online_query = toc(time_online_query);

period_total = toc(time_total);


fprintf('\n========== new_ppr_mut_inc ============\n');

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





