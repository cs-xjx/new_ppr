% Generate topological sort on pws under mutiple edge, for each pw node, find its children
% values differ at only one position between parent and child nodes
% let the position is p
% the value of each node differs from that of its parent only after p

% tar = {[0 3], [0 5 6 8], [0 4 8]};
target = {[1 2], [3 4 5], [6 7]};
% tar = {[1 2 34 5], [4 5 6 8 7 9], [7 4 8 5 3 1]};

tarNodes = cat(2, tar{:}); % tarNodes: nodes in target, except -1
tarNodes = tarNodes(tarNodes ~= 0);
ntar = size(tarNodes, 2);
npw = 2 ^ ntar;

res = zeros(npw, ntar); % res: the array that stores the results
inheritFrom = zeros(npw, 1); % the pw inherit from

p = 0; % p: the different index of the current node and its parent
parent = 1; % parent: point to a pw, and the alg generates the children of the current node
insertNew = 2; % insertNew: the index where the new node will be inserted

while parent < npw
    tempRes = res(parent, :);
    % only change values after p position as child nodes of the current node,
    % and only change the ith position every loop
    child_num = ntar - p;
    res(insertNew : insertNew + child_num - 1, :) = repmat(tempRes, child_num, 1);
    res(insertNew : insertNew + child_num - 1, p + 1 : ntar) = eye(child_num);
    inheritFrom(insertNew : insertNew + child_num - 1) = parent;
    insertNew = insertNew + child_num;
%     for i = p + 1 : ntar
%         tempRes(i) = 1 - tempRes(i);
%         res(insertNew, :) = tempRes;
%         inheritFrom(insertNew) = parent;
%         insertNew = insertNew + 1;
%         tempRes(i) = 1 - tempRes(i);
%     end
    parent = parent + 1;
    if parent >= insertNew
        break;
    end
    nextParentTar = res(parent, :);
    nextParentTarP = res(inheritFrom(parent), :);
    p = find(nextParentTar ~= nextParentTarP);
end
% res
% inheritFrom
