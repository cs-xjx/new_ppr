function T = flattenUnGraph_mul(A, source, target, n_certain)
%flattenUnGraph_mul  flatten each uncertain edge into certain edges(under multipel edge semantics)

new_nodes = cellfun(@(x) x(x > n_certain), target, 'UniformOutput', false);
new_nodes = cat(2, new_nodes{:});
n = size(A, 1); % n: # of nodes in graph

c = 0; % c: # of centain edges
u = 0; % u: # of uncentain edges
% globalRedistribute = 0;
signNode = sparse(1, n); % signNode: vector that marks the target node of the node
t = 2;

for i = 1 : n
    c = sum(A(i, :) ~= 0);
    ueIndexes = source == i; % index of the uncertain edges of a node
    tarNodes = cat(2, target{ueIndexes}); % tarNodes: nodes in target, except -1
    tarNodes = tarNodes(tarNodes ~= 0); % % change -1 into 0
    u = numel(tarNodes);
    if c == 0 && u == 0
        continue;
    end
    signNode(A(i, :) ~= 0) = 1;
    A(i, :) = A(i, :) / (c + u);
    u0 = 0; % # of uncertain edges which does not have epsilon(-1)
    ueTranPro = (1 / (c + u)) * (1 / t);
    for j = 1 : u
        ue = tarNodes(j);
        A(i, ue) = A(i, ue) + ueTranPro;
    end
    partRedistribute = ueTranPro * u;
    if partRedistribute == 0
        signNode(:) = 0;
        continue;
    end
    if c > 0
        signedNode = find(signNode == 1);
        A(i, signedNode) = A(i, signedNode) + partRedistribute / c;
        signNode(:) = 0;
    else
%         globalRedistribute = globalRedistribute + partRedistribute;
%         A(i, :) = A(i, :) + partRedistribute / n;
        A(i, [1 : n_certain new_nodes]) = A(i, [1 : n_certain new_nodes]) + partRedistribute / (n_certain + numel(new_nodes));
    end
end
% A(A ~= 0) = A(A ~= 0) + globalRedistribute / n;
T = A';

end

