function [res, inheritFrom, isleaf] = bfs_topologicalsort_mut(tar)
% Generate topological sort on pws under mutual exclusion, for each pw node, find its children
% values differ at only one position between parent and child nodes
% let the position is p
% the value of each node differs from that of its parent only after p

% tar = {[-1 3], [-1 5 6 8], [-1 4 8]};
% tar = {[4 7 2 8], [3 5 6 9], [6 9 4], [1 5], [9 2]};
% tar = {[1 2], [3 3 4 5], [6 7]};
% tar = {[4 7 2 8 5], [3 5 6 9 4], [6 9 4 3 2], [1 5 6 4 9], [9 2 6 7 3]};

ue_num = numel(tar); % numUE: # of uncertain edges
root = zeros(1, ue_num);
npw = 1;

% select the 1th value of each uncertain edge as the initial value of res
for i = 1 : ue_num
    root(i) = tar{i}(1);
    npw = npw * numel(tar{i});
end

res = zeros(npw, ue_num); % res: the array that stores the results
inheritFrom = zeros(npw, 1); % the pw inherit from
isleaf = ones(npw, 1);
diff = zeros(npw, 1);

res(1, :) = root;
% inheritFrom(1) = 0;
% isleaf(1) = 0;

p = 0; % p: the different index of the current node and its parent
parent = 1; % parent: point to a pw, and the alg generates the children of the current node
insertNew = 2; % insertNew: the index where the new node will be inserted

while parent < npw
    tempRes = res(parent, :);
    if  p < ue_num
        isleaf(parent) = 0;
    end
    % only change values after p position as child nodes of the current node,
    % and only change the ith position every loop
    for i = p + 1 : ue_num
%         dealTar = setdiff(target{i}, tempRes(i)); % dealTar: values can be replaced at the position
%         dealTar = tar{i}(tar{i} ~= tempRes(i));
        dealTar = tar{i};
        f = find(dealTar == tempRes(i));
        dealTar(f(1)) = [];
        deal_n = numel(dealTar);
        for j = 1 : deal_n
            tempRes(i) = dealTar(j);
            res(insertNew, :) = tempRes;
            inheritFrom(insertNew) = parent;
            diff(insertNew) = i;
            insertNew = insertNew + 1;
        end
        tempRes(i) = res(parent, i);
    end
    parent = parent + 1;
%     p = find(res(parent, :) ~= res(inheritFrom(parent), :));
    while parent < npw && diff(parent) == ue_num
        parent = parent + 1;
    end
    p = diff(parent);
end
% res;
% inheritFrom
