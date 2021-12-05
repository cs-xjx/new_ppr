function A = BearS_flat_mut(A, source, target, na)
%flattenUnGraph_mut  flatten each uncertain edge into certain edges(under mutual exclusion semantics)

n = size(A, 1); % n: # of nodes in graph
ue_num = size(source, 2); % numUE: # of uncertain edge

new_nodes = cellfun(@(x) x(x > na), target, 'UniformOutput', false);
new_nodes = cat(2, new_nodes{:});

unique_src = unique(source);
unique_n = numel(unique_src);

% globalRedistribute = 0;
signCerNode = sparse(1, n); % signNode: vector that marks the target node of the node
signU0Edge = sparse(1, ue_num);
for src_i = 1 : unique_n
    i = source(src_i);
    c = sum(A(i, :) ~= 0);
    ueIndexes = find(source == i); % index of the uncertain edges of a node
    u = numel(ueIndexes);

%     signNode(A(i, :) ~= 0) = 1;
    signCerNode = A(i, :);
    A(i, :) = A(i, :) / (c + u);
    u0 = 0; % # of uncertain edges which does not have epsilon(-1)
    partRedistribute = 0;
    for j = 1 : u
        ue = target{ueIndexes(j)};
        t = size(ue, 2);
        tranPro = (1 / (c + u)) * (1 / t);
        A(i, ue(ue ~= 0)) = A(i, ue(ue ~= 0)) + tranPro; % % change -1 into 0
        if all(ue ~= 0)
            u0 = u0 + 1;
            % signNode(ue) = 1;
            signU0Edge(ueIndexes(j)) = 1;
        else
            partRedistribute = partRedistribute + tranPro;
        end
    end
    if c > 0 || u0 > 0
        signedCerNode = find(signCerNode);
        A(i, signedCerNode) = A(i, signedCerNode) + partRedistribute / (c + u0);
        signCerNode(:) = 0;
        
        signedU0Edge = find(signU0Edge);
        for k = 1 : u0
            ue = target{signedU0Edge(k)};
            t = size(ue, 2);
            A(i, ue) = A(i, ue) + (partRedistribute / (c + u0)) * (1 / t);
        end
        signU0Edge(:) = 0;
    else
%         A(i, :) = A(i, :) + partRedistribute / n;
        A(i, [1 : na new_nodes]) = A(i, [1 : na new_nodes]) + partRedistribute / (na + numel(new_nodes));
    end
    minPro = min(A(i, (A(i, :) > 0)));
    A(i, :) = A(i, :) ./ minPro;
end

end

