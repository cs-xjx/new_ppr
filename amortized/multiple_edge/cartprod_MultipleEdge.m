function [pwG] = cartprod_MultipleEdge(target, ue_num)
%cartprod_MultipleEdge get cartesian product of uncertain edges under multiple edge semantics
%   target: target nodes of uncertain edge
%   ue_num: # of uncertain edge
%   numOfWorlds: # of all possible worlds(# of rows of pwG)

% target = {[1 2], [3 4 5], [6 7]};
% ue_num = numel(target);

% degreeOfUE = zeros(1, ue_num); % possOfUE: # of degrees of each uncertain edge
% for i = 1 : ue_num
%     degreeOfUE(i) = sum(target{i} ~= 0);
% end

% tic;
% degree = cellfun(@(x) numel(x(x ~= 0)), target);
% tarNodes = cat(2, target{:}); % tarNodes: nodes in target, except -1
% tarNodes = tarNodes(tarNodes ~= 0);
% exist = [0 1];
% numOfTarNodes = size(tarNodes, 2);
% 
% % % method 1, get cartprod by ndgrip()
% tempCell = cell(1, numOfTarNodes);
% % [tempCell{1 : numOfTarNodes}] = ndgrid(exist);
% [tempCell{end : -1 : 1}] = ndgrid(exist);
% pwG_01 = reshape(cat(numOfTarNodes, tempCell{:}), [], numOfTarNodes);
% npw = size(pwG_01, 1);
% pwG_01 = logical(pwG_01);
% pwG_01cell = mat2cell(pwG_01, ones(npw, 1), degree);
% pwG = cell(npw, ue_num);
% for i = 1 : ue_num
%     pwG(:, i) = cellfun(@(x) target{i}(x), pwG_01cell(:, i), 'UniformOutput', false);
% end
% toc
% 

% handle pwG
% if target{i} doesnt include epsilon, then the edge cannot be in a situation where all target nodes are not selected
% scope_start = 1;
% for i = 1 : ue_num
%     if any(target{i} == 0)
%         scope_start = scope_start + degreeOfUE(i);
%         continue;
%     end
%     check_all_zero = sum(pwG(:, scope_start : scope_start + degreeOfUE(i) - 1), 2);
%     pwG(check_all_zero == 0, :) = [];
%     scope_start = scope_start + degreeOfUE(i);
% end

% tic;
tempCell = cell(1, ue_num);
for i = 1 : ue_num
    tari = target{i};
    tari = tari(tari ~= 0);
    ni = numel(tari);
    npwi = 2 ^ ni;
    pwiCell = cell(npwi, 1);
    pwiCell{1} = [];
    index = 2;
    for k = 1 : ni
        cnk = nchoosek(tari, k);
        rows = size(cnk, 1);
        pwiCell(index : index + rows - 1) = mat2cell(cnk, ones(1, rows));
        index = index + rows;
    end
    tempCell{i} = pwiCell;
end
[tempCell{1 : ue_num}] = ndgrid(tempCell{:});
pwG = reshape(cat(ue_num, tempCell{:}), [], ue_num);
% toc;

end

