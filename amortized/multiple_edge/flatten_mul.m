function A = flatten_mul(A, src, tar, n_certain)
%flatten_mul  flatten each uncertain edge into certain edges(under multipel edge semantics)

new_nodes = cellfun(@(x) x(x > n_certain), tar, 'UniformOutput', false);
new_nodes = cat(2, new_nodes{:});

unique_src = unique(src);
unique_n = numel(unique_src);

t = 2;

for src_i = 1 : unique_n
    i = unique_src(src_i);
    c = sum(A(i, :) ~= 0);
    ueIndexes = src == i; % index of the uncertain edges of a node
    
    tar_nodes = cat(2, tar{ueIndexes}); % tarNodes: nodes in target, except -1
    tar_nodes = tar_nodes(tar_nodes ~= 0); % % change -1 into 0
    u = numel(tar_nodes);
    
    signCerNode = A(i, :);
    A(i, :) = A(i, :) / (c + u);
%     u0 = 0; % # of uncertain edges which does not have epsilon(0)
    ue_tran_pro = (1 / (c + u)) * (1 / t);
    for j = 1 : u
        ue = tar_nodes(j);
        A(i, ue) = A(i, ue) + ue_tran_pro;
    end
    redistribution = ue_tran_pro * u;
    if c > 0
        signedCerNode = find(signCerNode);
        A(i, signedCerNode) = A(i, signedCerNode) + redistribution / c;
    else
%         A(i, :) = A(i, :) + redistribution / n;
        A(i, [1 : n_certain new_nodes]) = A(i, [1 : n_certain new_nodes]) + redistribution / (n_certain + numel(new_nodes));
    end
end

end

