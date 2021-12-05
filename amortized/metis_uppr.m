function [r, temp_period_partAndInv] = metis_uppr(Tcertain, P, c, n_blocks, qu)
%metis_ppr

n = size(Tcertain, 1);
s = sparse(qu, 1, 1 / numel(qu), n, 1); % select nodes in qu as important nodes

% k = 3;  %the number of partitions
if n_blocks > n
    error('the value of k should not greater than n');
end
time_partAndInv = tic;
[part, ~] = metismex('PartGraphRecursive', sparse(Tcertain), n_blocks); % small n_blocks
% [part, ~] = metismex('PartGraphKway', sparse(Tcertain), k); % large n_blocks

[~, nodeOrder] = sort(part); % reorder after partitioning£¬record sort order
% kNum = hist(part, unique(part)); % record the number of nodes for each partition
orderDiff = diff(nodeOrder);
ltzeroIndex = find(orderDiff < 0);
kNum = diff(ltzeroIndex);
kNum = [ltzeroIndex(1) kNum n - ltzeroIndex(end)];

% get the reordered matrix, reorderTc
reorderTc = Tcertain(nodeOrder, nodeOrder);

% get TBL, invQBL and TX
k2 = 0;
TBLk = cell(n_blocks, 1);
invQBLk = cell(n_blocks, 1);
for i = 1 : n_blocks
    k1 = k2 + 1;
    k2 = k2 + kNum(i);
    TBLk{i} = reorderTc(k1 : k2, k1 : k2);
    invQBLk{i} = inv(speye(kNum(i)) - c * TBLk{i});
%     TBLk = reorderTc(k1 : k2, k1 : k2);
%     invQBLk = inv(speye(kNum(i)) - c * TBLk);
%     if i == 1
%         TBL = TBLk;
%         invQBL = invQBLk;
%     else
%         TBL = blkdiag(TBL, TBLk);
%         invQBL = blkdiag(invQBL, invQBLk);
%     end
end
TBL = blkdiag(TBLk{:});
invQBL = blkdiag(invQBLk{:});
TX = reorderTc - TBL;
temp_period_partAndInv = toc(time_partAndInv);

P = P(nodeOrder, nodeOrder);
s = s(nodeOrder);

qs = invQBL * s; %1
tqs = TX * qs; %2
qtqs = invQBL * tqs; %2
pqs = P * qs; %3
qpqs = invQBL * pqs; %3
tqtqs = TX * qtqs; %4
qtqtqs = invQBL * tqtqs; %4
pqtqs = P * qtqs; %5
qpqtqs = invQBL * pqtqs; %5
tqpqs = TX * qpqs; %6
qtqpqs = invQBL * tqpqs; %6

r = (1-c)*qs + c*(1-c)*qtqs + (c*(1-c))*qpqs + c*c*(1-c)*qtqtqs + (c*c*(1-c))*qpqtqs + (c*c*(1-c))*qtqpqs;

% r = (1 - c) * (speye(n) + (c / pw_n) * invQBL * ((pw_n * TX + P) + ...
%                c * (pw_n * TX * invQBL * TX + P * invQBL * TX + TX * invQBL * P))) * invQBL * s;
[~, nodeOrder2] = sort(nodeOrder);
r = r(nodeOrder2);

end

