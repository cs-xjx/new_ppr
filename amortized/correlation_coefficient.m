function mrho = correlation_coefficient(standard_pprs, target_pprs, type)
%correlation_coefficient 此处显示有关此函数的摘要
% 返回平均相关系数

mrho = 0;

qu_num = numel(standard_pprs);
for qi = 1 : qu_num
    rho = corr(standard_pprs{qi}, target_pprs{qi}, 'Type', type);
    mrho = mrho + rho;
end
mrho = mrho / qi;

end

