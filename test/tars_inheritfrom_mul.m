function [src_index, tar_nodes, inheritFrom, level, diff_p, isleaf] = tars_inheritfrom_mul(src, tar)
%tars_inheritfrom_mul  get cartesian product and inherit relationship of uncertain edges under multiple edge semantics
%   tar: target nodes of uncertain edge
%   ue_num: # of uncertain edge
%   numOfWorlds: # of all possible worlds(# of rows of pwG)

% src = [4 6 4];
% tar = {[1 2], [3 4 5], [6 7]};
% ue_num = numel(target);

[sorted_src, idx] = sort(src);
tar = tar(idx);

degree = cellfun('length', tar);
tar_nodes = cat(2, tar{:}); % tarNodes: nodes in target, except -1
ntar = size(tar_nodes, 2);
src_index = ones(1, ntar);

start = degree(1) + 1;
for i = 2 : numel(sorted_src)
    if sorted_src(i) ~= sorted_src(i - 1)
        src_index(start : end) = src_index(start : end) + 1;
    end
    start = start + degree(i);
end

npw = 2 ^ ntar;


parent_index = zeros(ntar, 1);
parent_index(1) = 1;
pos = 1;

diff_p = zeros(npw, 1);
diff_c = ones(npw, 1) * ntar;

inheritFrom = zeros(npw, 1);
level = zeros(npw, 1);
isleaf = zeros(npw, 1);
level(1) = 1;

for i = 2 : npw
    parent = parent_index(pos);
    inheritFrom(i) = parent;
    diff_p(i) = diff_c(parent);
    level(i) = level(parent) + 1;
    
    diff_c(parent) = diff_c(parent) - 1;
    if diff_c(parent) == diff_p(parent)
        parent_index(pos) = 0;
    end
    
    if diff_p(i) ~= ntar
        pos = pos + 1;
        parent_index(pos) = i;
    else
        isleaf(i) = 1;
        while pos > 0 && parent_index(pos) == 0
            pos = pos - 1;
        end
    end
end




