function [res, inheritFrom, isleaf, level, diff_p] = dfs_topologicalsort_mut(tar)
% Generate topological sort on pws under mutual exclusion, for each pw node, find its children
% values differ at only one position between parent and child nodes
% let the position is p
% the value of each node differs from that of its parent only after p

% target = {[-1 3], [-1 5 6 8], [-1 4 8]};
% tar = {[4 7 2 8], [3 5 6 9], [6 9 4], [1 5], [9 2]};
% tar = {[1 2], [3 4 5], [6 7]};

ue_num = numel(tar); % numUE: # of uncertain edges
root = zeros(1, ue_num);
npw = 1;

% select the 1th value of each uncertain edge as the initial value of res
for i = 1 : ue_num
    root(i) = tar{i}(1);
    npw = npw * numel(tar{i});
end

res = zeros(npw, ue_num); % res: the array that stores the results
inheritFrom = zeros(npw, 1);
isleaf = ones(npw, 1);
level = zeros(npw, 1);
diff_p = zeros(npw, 1);

res(1, :) = root;
inheritFrom(1) = 0;
isleaf(1) = 0;
level(1) = 1;

% p = 0; % p: the different index of the current node and its parent
% parent = 1; % parent: point to a pw, and the alg generates the children of the current node
insertNew = 2; % insertNew: the index where the new node will be inserted

temp_root = root;
for i = 1 : ue_num
    deal_tar = tar{i};
    f = find(deal_tar == root(i));
    deal_tar(f(1)) = [];
    deal_n = numel(deal_tar);
    for j = 1 : deal_n
        temp_root(i) = deal_tar(j);
        res(insertNew, :) = temp_root;
        inheritFrom(insertNew) = 1;
        level(insertNew) = 2;
        diff_p(insertNew) = i;
        
        insertNew = insertNew + 1;
        dfs(insertNew - 1, i);
    end
    temp_root(i) = root(i);
end

% res
% inheritFrom
% level

function dfs(parent, diff)
    temp_parent_tar = res(parent, :);
    if diff < ue_num
        isleaf(parent) = 0;
    end
    for p = diff + 1 : ue_num
        deal_tar2 = tar{p};
        f2 = find(deal_tar2 == temp_parent_tar(p));
        deal_tar2(f2(1)) = [];
        deal_n2 = numel(deal_tar2);
        for k = 1 : deal_n2
            temp_parent_tar(p) = deal_tar2(k);
            res(insertNew, :) = temp_parent_tar;
            inheritFrom(insertNew) = parent;
            level(insertNew) = level(parent) + 1;
            diff_p(insertNew) = p;
        
            insertNew = insertNew + 1;
            dfs(insertNew - 1, p);
        end
        temp_parent_tar(p) = res(parent, p);
    end
end


% while parent < insertNew
%     tempRes = res(parent, 1 : ue_num);
%     if  p < ue_num
%         res(parent, ue_num + 2) = 0;
%     end
%     % only change values after p position as child nodes of the current node,
%     % and only change the ith position every loop
%     for i = p + 1 : ue_num
%         dealTar = tar{i};
%         f = find(dealTar == tempRes(i));
%         dealTar(f(1)) = [];
%         deal_n = numel(dealTar);
%         for j = 1 : deal_n
%             tempRes(i) = dealTar(j);
%             res(insertNew, :) = [tempRes parent 1];
%             insertNew = insertNew + 1;
%         end
%         tempRes(i) = res(parent, i);
%     end
%     parent = parent + 1;
%     if parent >= insertNew
%         break;
%     end
%     p = find(res(parent, 1 : ue_num) ~= res(res(parent, ue_num + 1), 1 : ue_num));
% end

end
