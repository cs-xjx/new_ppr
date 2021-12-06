function [inheritFrom, level, diff_p, isleaf] = inheritance_relation(tar, npw)
%      get inheritance relationship of tars

nsrc = numel(tar);
child_num = cellfun(@(x) numel(x) - 1, tar);
parent_index = zeros(nsrc, 1);
parent_index(1) = 1;
pos = 1;

diff_p = zeros(npw, 1);
diff_c = ones(npw, 1) * nsrc;
has_child_now = zeros(npw, 1);

inheritFrom = zeros(npw, 1);
level = zeros(npw, 1);
isleaf = zeros(npw, 1);
level(1) = 1;

for i = 2 : npw
    parent = parent_index(pos);
    inheritFrom(i) = parent;
    diff_p(i) = diff_c(parent);
    level(i) = level(parent) + 1;
    
    now_child = has_child_now(parent) + 1;
    if now_child < child_num(diff_c(parent))
        has_child_now(parent) = now_child;
    else
        diff_c(parent) = diff_c(parent) - 1;
        has_child_now(parent) = 0;
        if diff_c(parent) == diff_p(parent)
            parent_index(pos) = 0;
        end
    end
    if diff_p(i) ~= nsrc
        pos = pos + 1;
        parent_index(pos) = i;
    else
        isleaf(i) = 1;
        while pos > 0 && parent_index(pos) == 0
            pos = pos - 1;
        end
    end
end

end
