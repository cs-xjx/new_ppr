% function main_mul()

clear all;
addpath('D:\Polyspace\metis-5.1.0\metismex-dgleich\');
% addpath('D:\Polyspace\metis-5.1.0\metismex-dgleich\yu\');
addpath('D:\Polyspace\metis-5.1.0\metismex-dgleich\amortized\multiple_edge\');
addpath('D:\Polyspace\metis-5.1.0\metismex-dgleich\test\');
% addpath('D:\Polyspace\metis-5.1.0\metismex-dgleich\BEAR-2.0\');

% get a
fpath = 'D:\Polyspace\metis-5.1.0\metismex-dgleich\dataset\';
%fpath = 'E:\datasets\Matlab\';
% ds = '494_bus';
% ds = 'ego-Facebook';
ds = 'wiki-Vote';
% ds = 'p2p-Gnutella31';
% ds = 'email-EuAll';
% ds = 'web-NotreDame';
% ds = 'dblp_03-05';
fname = [fpath, ds, '.mat'];
load(fname);
a = Problem.A;
clear Problem;
%
% datasets: https://www.cise.ufl.edu/research/sparse/matrices/SNAP/
%
% as-735 dblp_03-05 dblp_06-08 dblp_09-11
% cit-HepPh p2p-Gnutella08 email-Enron cit-Patents web-Stanford 
% web-BerkStan web-Google email-EuAll uk-2002 soc-LiveJournal1


% get src, tar, qu
src = [351 5147 3348 1010];
tar = {[1576 4871 3001], [160 8414 1920], [8330 2745 8339], [783 3313 8339]};
qu = {[5147 4196 6129 6315 2100 2359 1505 1010 351 6560], ...
    [5147 7384 6604 7915 1213 1641 1403 1010 3348 8], ...
    [5147 603 4826 6660 2776 1626 4705 1010 351 87], ...
    [5147 1265 3400 610 4382 3827 8189 1010 3348 1972], ...
    [5147 4293 7866 7801 6188 4925 4737 1010 351 7231]};


% set c, kmax, num_blks, low rank
c = 0.6;            % decay factor
kmax = 20;         % # of iterations
n_blocks = 50;        % # of partitions
r = 20;        % low rank of svds


%qu_len = 100;       % query length
%qu_num = 10;        % # of queries
%ue_num = 4;  % 8;   % # of uncertain edges (source nodes)
%ue_deg = 10; % 4;   % # of uncertain target nodes for each uncertain edge

uede = ['-', num2str(numel(src)), '_', num2str(numel(tar{1}) + 1)];
% algo
% [res1] = exhPPR_MultipleEdge(a, c, kmax, src, tar, qu, ds);
% save(['exh_mul-', ds, uede], 'res1');

standard_res  = load(['exh_mul-', ds, uede]);
res1 = standard_res.res1;

% [res2] = exhPPR(a, c, kmax, src, tar, qu, ds);

% [res_ea] = exhApxPPR_MultipleEdge(a, c, n_blocks, src, tar, qu, r, ds);
% save(['eap_mul-', num2str(n_blocks), ds, uede], 'res_ea');

[res_cp] = collPPR_MultipleEdge(a, c, kmax, src, tar, qu, ds);
save(['cp_mul-', ds, uede], 'res_cp');

% [res_cap] = collApxPPR_MultipleEdge(a, c, n_blocks, src, tar, qu, r, ds);
% save(['cap_mul-', num2str(n_blocks), ds, uede], 'res_cap');

% [res_ca2p] = collApx2PPR_MultipleEdge(a, c, n_blocks, src, tar, qu, r, ds);
% save(['ca2p_mul-', ds, uede], 'res_ca2p');

% [res_fp] = flatPPR_MultipleEdge(a, c, kmax, src, tar, qu, ds);
% save(['fp_mul-', ds, uede], 'res_fp');

% [res_fap] = flatApxPPR_MultipleEdge(a, c, n_blocks, src, tar, qu, r, ds);
% save(['fap_mul-', num2str(n_blocks), ds, uede], 'res_fap');

% [res_bear] = BearS_uncertain_mul(a, c, src, tar, qu, ds);
% save(['bear_mul-', ds, uede], 'res_bear');

% [res_uppr] = UPPR_MultipleEdge(a, c, src, tar, qu, n_blocks,ds);
% save(['uppr_mul-', num2str(n_blocks), ds, uede], 'res_uppr');

% [res3] = new_ppr_mul_inc(a, c, kmax, src, tar, qu, ds);


% check accuracy
top_k = 50;
[sum_err, mAP, nDCG, kendall, spearman] = compute_precision(res1, res_cp, top_k, 1);
fprintf('\n\n========== accuracy ============\n');
sum_err
mAP
nDCG
kendall
spearman
fprintf(' >>         error              :  %f   \n', sum_err);
fprintf(' >>         mAP                :  %f   \n', mAP);
fprintf(' >>         nDCG               :  %f   \n', nDCG);
fprintf(' >>         kendall            :  %f   \n', kendall);
fprintf(' >>         spearman           :  %f   \n', spearman);


% check accuracy
% qu_num = numel(qu);
% sum_err = 0;
% for qi = 1:qu_num
%     sum_err = sum_err + norm(res1{qi} - res_cap{qi});
% end
% fprintf('\n\n========== accuracy ============\n');
% sum_err
% fprintf(' >>         error              :  %f  \n\n\n', sum_err);



