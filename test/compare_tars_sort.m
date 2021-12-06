% tar = {[3536 833 2593], [4324 4071 4089], [210 2948 256], [3977 3149 712], [4049 1209 1337], [102 2260 1405], [4093 1336 477], [2854 1162 1877], [404 4147 4196], [4102 2733 3686]};
tar = {[1 2], [3 4 5], [6 7]};
% tar = {[4 7 2], [3 5 6 9], [6 9 4], [1 5], [9 2]};
nsrc = numel(tar);
tar0 = cellfun(@(x) [0 x], tar, 'UniformOutput', false);
tar0 = tar;

% 
% 
% 0. use previous sort
% 
fprintf(' >> start method 0, use previous sort ... \n');
time_gen_pwd = tic;
% tar0 = cellfun(@(x) [0 x], tar, 'UniformOutput', false);
% tar0 = tar;
tars = cell(1,nsrc);
[tars{end:-1:1}] = ndgrid(tar0{end:-1:1});
tars = cat(nsrc + 1, tars{:});
tars = reshape(tars, [], nsrc);  % tars: possible world
npw = size(tars, 1);
period_gen_pwd = toc(time_gen_pwd);
fprintf('  # of possible worlds            :  %d  \n', npw);
fprintf('  Time (generate possible worlds) :  %fs \n\n', period_gen_pwd);
%
%
clear tars npw;

% 
% 
% 1. use bfs topological sort
% 
fprintf(' >> start method 1, use bfs topological sort ... \n');
time_gen_bfs_pwd = tic;
% tar0 = tar;
[tars_bfs, inheritFrom, isleaf] = bfs_topologicalsort_mut(tar0);
npw = size(tars_bfs, 1);
period_gen_bfs_pwd = toc(time_gen_bfs_pwd);
fprintf('  # of possible worlds            :  %d  \n', npw);
fprintf('  Time (generate bfs sort) :  %fs \n\n', period_gen_bfs_pwd);
%
%
clear tars_bfs npw inheritFrom isleaf;

% 
% 
% 2. use dfs topological sort
% 
fprintf(' >> start method 2, use dfs topological sort ... \n');
time_gen_dfs_pwd = tic;
% tar0 = tar;
[tars_dfs, inheritFrom, isleaf, level] = dfs_topologicalsort_mut(tar0);
npw = size(tars_dfs, 1);
period_gen_dfs_pwd = toc(time_gen_dfs_pwd);
fprintf('  # of possible worlds            :  %d  \n', npw);
fprintf('  Time (generate dfs sort) :  %fs \n\n', period_gen_dfs_pwd);
%
%
clear tars_dfs npw inheritFrom isleaf level;

% 
% 
% 3. optimize dfs sort
% 
fprintf(' >> start method 3, optimize dfs sort ... \n');
time_gen_pwd = tic;
% tar0 = tar;
tars = cell(1,nsrc);
[tars{end:-1:1}] = ndgrid(tar0{end:-1:1});
tars = cat(nsrc + 1, tars{:});
tars = reshape(tars, [], nsrc);  % tars: possible world
npw = size(tars, 1);
% 
% child_num = zeros(nsrc, 1);
% for i = 1 : nsrc
%     child_num(i) = numel(tar{i}) - 1;
% end
child_num = cellfun(@(x) numel(x) - 1, tar0);
parent_index = zeros(nsrc, 1);
parent_index(1) = 1;
pos = 1;

diff_p = zeros(npw, 1);
diff_c = ones(npw, 1) * nsrc;
has_child_now = zeros(npw, 1);

inheritFrom = zeros(npw, 1);
level = zeros(npw, 1);
isleaf = zeros(npw, 1);
level(1) = 1;

for i = 2 : npw
    parent = parent_index(pos);
    inheritFrom(i) = parent;
    diff_p(i) = diff_c(parent);
    level(i) = level(parent) + 1;
    
    now_child = has_child_now(parent) + 1;
    if now_child < child_num(diff_c(parent))
        has_child_now(parent) = now_child;
    else
        diff_c(parent) = diff_c(parent) - 1;
        has_child_now(parent) = 0;
        if diff_c(parent) == diff_p(parent)
            parent_index(pos) = 0;
        end
    end
    if diff_p(i) ~= nsrc
        pos = pos + 1;
        parent_index(pos) = i;
    else
        isleaf(i) = 1;
        while pos > 0 && parent_index(pos) == 0
            pos = pos - 1;
        end
    end
%     if pos <= 0 && i < npw
%         i
%     end
end
% 
period_gen_pwd = toc(time_gen_pwd);
fprintf('  # of possible worlds            :  %d  \n', npw);
fprintf('  Time (generate possible worlds) :  %fs \n\n', period_gen_pwd);
%
%

% % % check if the vector inheritFrom of method 3 we get is true
istrue = 'true';
for i = 2 : npw
    parent = inheritFrom(i);
    if parent == i
        istrue = 'false1';
        break;
    end
    diff = find(tars(parent, :) ~= tars(i, :));
    if numel(diff) > 1 || (numel(diff) == 1 && diff ~= diff_p(i))
        istrue = 'false2';
        break;
    end
    if (diff == nsrc && isleaf(i) == 0) || (diff ~= nsrc && isleaf(i) == 1)
        istrue = 'false3';
        break;
    end
    if parent < i - 1 && parent > inheritFrom(i - 1)
        istrue = 'false4';
        break;
    end
    if level(i) ~= level(parent) + 1
        istrue = 'false5';
        break;
    end
end
fprintf('  inheritFrom of method 3     :  %s  \n\n', istrue);

clear all;





