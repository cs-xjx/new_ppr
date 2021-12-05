function A = flatten_mut(A, src, tar, n_certain)
%flatten_mut  flatten each uncertain edge into certain edges(under mutual exclusion semantics)

new_nodes = cellfun(@(x) x(x > n_certain), tar, 'UniformOutput', false);
new_nodes = cat(2, new_nodes{:});

unique_src = unique(src);
unique_n = numel(unique_src);

for src_i = 1 : unique_n
    i = unique_src(src_i);
    c = sum(A(i, :) ~= 0);
    ueIndexes = find(src == i); % index of the uncertain edges of a node
    u = numel(ueIndexes);

%     signNode(A(i, :) ~= 0) = 1;
    signCerNode = A(i, :);
    A(i, :) = A(i, :) / (c + u);
    u0 = 0; % # of uncertain edges which does not have epsilon(0)
    partRedistribute = 0;
    for j = 1 : u
        ue = tar{ueIndexes(j)};
        t = size(ue, 2);
        tranPro = (1 / (c + u)) * (1 / t);
        A(i, ue(ue ~= 0)) = A(i, ue(ue ~= 0)) + tranPro; % % change -1 into 0

        partRedistribute = partRedistribute + tranPro;
    end
    if c > 0 || u0 > 0
        signedCerNode = find(signCerNode);
        A(i, signedCerNode) = A(i, signedCerNode) + partRedistribute / (c + u0);
        
    else
%         A(i, :) = A(i, :) + partRedistribute / n;
        A(i, [1 : n_certain new_nodes]) = A(i, [1 : n_certain new_nodes]) + partRedistribute / (n_certain + numel(new_nodes));
    end
%     minPro = min(A(i, A(i, :) > 0));
%     A(i, :) = A(i, :) ./ minPro;
end

end

