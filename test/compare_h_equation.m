% function [res] = new_ppr_mut(a, c, kmax, src, tar, qu, ds)
clear all;

% get a
fpath = 'D:\Polyspace\metis-5.1.0\metismex-dgleich\dataset\';
%fpath = 'E:\datasets\Matlab\';
% ds = '494_bus';
ds = 'ego-Facebook';
% ds = 'p2p-Gnutella31';
% ds = 'email-EuAll';
% ds = 'web-NotreDame';
% ds = 'dblp_03-05';
fname = [fpath, ds, '.mat'];
load(fname);
a = Problem.A;
clear Problem;
src = [1981 1685 2509 3476 1610 1887 1802 3993 3086 1357];
tar = {[3536 833 2593], [4324 4071 4089], [210 2948 256], [3977 3149 712], [4049 1209 1337], [102 2260 1405], [4093 1336 477], [2854 1162 1877], [404 4147 4196], [4102 2733 3686]};
qu = {[260 553], [2658 145], [1163 3525], [2930 3248], [1457 1100]};

% src = [1059 3586 3572 3437];
% tar = {[4311 4199 1452 1598 2847 79 1719 795 1743], [3104 3263 1523 4324 1321 1771 3102 4350 2076], [623 4186 3323 1283 362 4089 1994 4218 4059], [3755 2353 3550 3 2474 2130 3236 4219 2319]};
% qu = {[260 553], [2658 145], [1163 3525], [2930 3248], [1457 1100]};

a = sparse([
0, 0, 0, 0, 1, 0, 0, 0, 0;
0, 0, 1, 1, 0, 0, 0, 0, 0;
0, 0, 0, 0, 0, 1, 0, 0, 0;
0, 0, 0, 0, 1, 0, 0, 0, 1;
0, 1, 0, 0, 0, 0, 1, 0, 0;
0, 1, 0, 0, 0, 0, 0, 0, 0;
0, 0, 1, 0, 0, 0, 0, 0, 0;
0, 0, 0, 1, 0, 0, 0, 0, 0;
1, 0, 0, 0, 0, 0, 0, 0, 0
    ]);
src = [1 2 1 2 5];
tar = {[4 7 2 8], [3 5 6 9], [6 9 4], [1 5], [9 2]};
qu = {[1 2 6]};

% time_total = tic;
fprintf('\n >> Start new_ppr_mut ... \n\n');

n_certain = size(a, 1); % n_certain: # of nodes in graph
% m_certain = nnz(a);

qu_num = numel(qu);
nsrc = numel(src);

n = max([n_certain, tar{:}]);

a(n, n) = 0;
%%%%%%%%%
c = 0.8;
%%%%%%%%%



suma = sum(a, 2);
d = 1 ./ suma;
d(~isfinite(d)) = 0;
Q = a' * spdiags(d, 0, n, n); % transition matrix of a
% D = spdiags(suma(src), 0, nsrc, nsrc);    % outdegree of source nodes
D = suma(src);
Dii = diag(D);

R = inv(eye(n) - c * Q);
En = speye(n);
El = speye(nsrc);

unique_src = unique(src);
[sorted_src, node_idx] = sort(src);
un_src_first = find(diff([-1 sorted_src n + 1]));
src_size = diff(un_src_first);
% src_index_map = containers.Map('KeyType','double','ValueType','any');
src_index_cell = cell(numel(unique_src), 1);
for i = 1 : numel(unique_src)
    scope = [un_src_first(i) : un_src_first(i + 1) - 1];
%     src_index_map(unique_src(i)) = node_idx(scope);
    src_index_cell{unique_src(i)} = node_idx(scope);
end

C = speye(nsrc);
dL = zeros(nsrc, 1);
for i = 1 : numel(unique_src)
    src_i = unique_src(i);
%     dL(src_index_map(src_i)) = src_size(i);
    dL(src_index_cell{src_i}) = src_size(i);
end
L = spdiags(dL, 0, nsrc, nsrc);


% 
% 
% 0. use previous sort
% 
fprintf(' >> start method 0, use previous sort ... \n');
%
% generate possible world (mutual exclusive semantics)
%
time_gen_pwd = tic;
fprintf(' >> Generate Possible World ... \n');
% tar0 = cellfun(@(x) [0 x], tar, 'UniformOutput', false);
tar0 = tar;
tars = cell(1,nsrc);
[tars{end:-1:1}] = ndgrid(tar0{end:-1:1});
tars = cat(nsrc + 1, tars{:});
tars = reshape(tars, [], nsrc);  % tars: possible world
npw = size(tars, 1);
period_gen_pwd = toc(time_gen_pwd);
fprintf('  # of possible worlds            :  %d  \n', npw);
fprintf('  Time (generate possible worlds) :  %fs \n', period_gen_pwd);
%
%
time_gen_h_0 = tic;
sum_h0 = zeros(nsrc, 1);
sum_H0 = zeros(nsrc, nsrc);
for qi = 1 : qu_num
    s = sparse(qu{qi}, 1, 1 / numel(qu{qi}), n, 1);
    p0 = R * s;
    b = C * p0(src);
    
    for i = 1 : npw
        tar_i = tars(i, :);
        H = inv(Dii + L - C * (c * R(src, tar_i) - R(src, src) + En(src, src)));
        h = H * b;
        sum_H0 = sum_H0 + H;
        sum_h0 = sum_h0 + h;
    end
end
period_gen_h_0 = toc(time_gen_h_0);
fprintf('  Time (generate h) :  %fs \n\n', period_gen_h_0);

% 
% 
% 1. use bfs topological sort
% 
fprintf(' >> start method 1, use bfs topological sort ... \n');

time_gen_bfs_pwd = tic;
[tars, inheritFrom, isleaf] = bfs_topologicalsort_mut(tar0);
npw = size(tars, 1);
period_gen_bfs_pwd = toc(time_gen_bfs_pwd);
fprintf('  # of possible worlds            :  %d  \n', npw);
fprintf('  Time (generate bfs sort) :  %fs \n', period_gen_bfs_pwd);

time_gen_h_1 = tic;
% H_map = containers.Map('KeyType','double','ValueType','any');
H_cell = cell(npw, 1);
h_cell = cell(npw, 1);
sum_h1 = zeros(nsrc, 1);
sum_H1 = zeros(nsrc, nsrc);
for qi = 1 : qu_num
    s = sparse(qu{qi}, 1, 1 / numel(qu{qi}), n, 1);
    p0 = R * s;
    b = C * p0(src);
    
    tar_i = tars(1, :);
    B = Dii + L - C * En(src, src) + C * R(src, src);
    H = inv(B - c * C * R(src, tar_i));
    h = H * b;
%     H_map(1) = H;
    H_cell{1} = H;
    h_cell{1} = h;
%     clear H;
    sum_H1 = sum_H1 + H;
    sum_h1 = sum_h1 + h;
    
    for i = 2 : npw
        tar_i = tars(i, :);
        parent = inheritFrom(i);
        parent_tar_i = tars(parent, :);
        j = find(parent_tar_i ~= tar_i);
        u = R(src, tar_i(j)) - R(src, parent_tar_i(j));
%         H1 = H_map(parent);
        H1 = H_cell{parent};
        y = H1 * (C * u);
        
        if isleaf(i)
%             h1 = h_map(parent);
            H2 = H1 + c / (1 - c * y(j)) * y * H1(j, :);
            h1 = H1 * b;
%             h1 = h_cell{parent};
            h = h1 + (c * h1(j)) / (1 - c * y(j)) * y;
%             ch1 = c * h1(src_j);
%             cy = 1 - c * y(src_j);
%             k = ch1 / cy;
%             delta_h = k * y;
%             h = h1 + delta_h;
        else
            H2 = H1 + c / (1 - c * y(j)) * y * H1(j, :);
            h = H2 * b;
%             H_map(i) = H2;
            H_cell{i} = H2;
%             h_cell{i} = h;
%             h_map(i) = h;
        end
        sum_H1 = sum_H1 + H2;
        sum_h1 = sum_h1 + h;
        
    end
end
period_gen_h_1 = toc(time_gen_h_1);
fprintf('  Time (generate h) :  %fs \n\n', period_gen_h_1);


% 
% 
% 2. use dfs topological sort
% 
fprintf(' >> start method 2, use dfs topological sort ... \n');

time_gen_dfs_pwd = tic;
[tars, inheritFrom, isleaf, level] = dfs_topologicalsort_mut(tar0);
npw = size(tars, 1);
period_gen_dfs_pwd = toc(time_gen_dfs_pwd);
fprintf('  # of possible worlds            :  %d  \n', npw);
fprintf('  Time (generate dfs sort) :  %fs \n', period_gen_dfs_pwd);

time_gen_h_2 = tic;
H_cell = cell(nsrc, 1);
h_cell = cell(nsrc, 1);
sum_h2 = zeros(nsrc, 1);
sum_H2 = zeros(nsrc, nsrc);
for qi = 1 : qu_num
    s = sparse(qu{qi}, 1, 1 / numel(qu{qi}), n, 1);
    p0 = R * s;
    b = C * p0(src);
    
    tar_i = tars(1, :);
    B = Dii + L - C * En(src, src) + C * R(src, src);
    H = inv(B - c * C * R(src, tar_i));
    h = H * b;
%     H_map(1) = H;
    H_cell{1} = H;
%     h_cell{1} = h;
%     clear H;
    sum_h2 = sum_h2 + h;
    sum_H2 = sum_H2 + H;
    
    for i = 2 : npw
        tar_i = tars(i, :);
        parent = inheritFrom(i);
        parent_tar_i = tars(parent, :);
        j = find(parent_tar_i ~= tar_i);
        u = R(src, tar_i(j)) - R(src, parent_tar_i(j));
        H1 = H_cell{level(parent)};
        y = H1 * (C * u);
        
        if isleaf(i)
            H2 = H1 + c / (1 - c * y(j)) * y * H1(j, :);
            h1 = H1 * b;
%             h1 = h_cell{level(parent)};
            h = h1 + (c * h1(j)) / (1 - c * y(j)) * y;
        else
            H2 = H1 + c / (1 - c * y(j)) * y * H1(j, :);
            h = H2 * b;
            H_cell{level(i)} = H2;
%             h_cell{level(i)} = h;
        end
        sum_h2 = sum_h2 + h;
        sum_H2 = sum_H2 + H2;
        
    end
end
period_gen_h_2 = toc(time_gen_h_2);
fprintf('  Time (generate h) :  %fs \n\n', period_gen_h_2);


sum_err1 = norm(sum_h0 - sum_h1);
sum_err2 = norm(sum_h0 - sum_h2);
fprintf('\n\n========== accuracy ============\n');
fprintf(' >>         error              :  %f  %f  \n\n\n', sum_err1, sum_err2);

