function mAP = mean_Average_Precision(standard_pprs, target_pprs, top_k)
%mean_Average_Precision 求向量的mAP

if (numel(standard_pprs) ~= numel(target_pprs))
    error('向量数目不相同');
end
% if top_k < 1
%     error('top_k needs > 0');
% end
all = 0;
if top_k == -1
    all = 1;
end
mAP = 0;
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
    precision = double(tag);
    for i = 2 : top_k
        precision(i) = precision(i) + precision(i - 1);
    end
    precision = precision ./ (1 : top_k)';
    max_precision = precision;
    for i = top_k - 1 : 1
        if max_precision(i) < max_precision(i + 1)
            max_precision(i) = max_precision(i + 1);
        end
    end
    AP = sum(max_precision .* tag) / top_k;
    mAP = mAP + AP;
end
mAP = mAP / qi;

end

