function [sum_err, mAP, nDCG, kendall, spearman] = compute_precision(res1, res2, top_k, compute_corr)

mAP = mean_Average_Precision(res1, res2, top_k);
nDCG = NDCG(res1, res2, top_k);
kendall = -1;
spearman = -1;
if compute_corr == 1
    kendall = correlation_coefficient(res1, res2, 'Kendall');
    spearman = correlation_coefficient(res1, res2, 'Spearman');
end



qu_num = numel(res1);
sum_err = 0;
for qi = 1 : qu_num
    sum_err = sum_err + norm(res1{qi} - res2{qi});
end
% fprintf('\n\n========== accuracy ============\n');
% sum_err
% mAP
% nDCG
% kendall
% spearman
% fprintf(' >>         error              :  %f   \n\n', sum_err);
% fprintf(' >>         mAP                :  %f   \n', mAP);
% fprintf(' >>         nDCG               :  %f   \n', nDCG);
% fprintf(' >>         kendall            :  %f   \n', kendall);
% fprintf(' >>         spearman           :  %f   \n', spearman);


