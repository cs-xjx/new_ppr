function nDCG = NDCG(standard_pprs, target_pprs, top_k)
%NDCG 此处显示有关此函数的摘要

if (numel(standard_pprs) ~= numel(target_pprs))
    error('向量数目不相同');
end
all = 0;
if top_k == -1
    all = 1;
end
nDCG = 0;
qu_num = numel(standard_pprs);
for qi = 1 : qu_num
    sppr = standard_pprs{qi};
    tppr = target_pprs{qi};
    if (numel(sppr) ~= numel(tppr))
        error(['向量', num2str(qi), '元素数目不相同']);
    end
    if top_k > numel(sppr)
        error(['top_k too large, ppr', num2str(qi), '''s num is ', num2str(numel(sppr))]);
    end
    if all == 1
        top_k = numel(sppr);
    end
    [~, s_idx] = sort(sppr, 'descend');
    [~, t_idx] = sort(tppr, 'descend');
    s_idx_topk = s_idx(1 : top_k);
    t_idx_topk = t_idx(1 : top_k);
    
    tag = s_idx_topk == t_idx_topk;
    tag = double(tag);
    s_ones = ones(top_k, 1);
    
    idcg = DCG(s_ones);
    dcg = DCG(tag);
    
    nDCG = nDCG + dcg / idcg;
end
nDCG = nDCG / qi;

    function dcg = DCG(s)
        dcg = sum(s ./ log2((1 : numel(s))' + 1));
    end

end

