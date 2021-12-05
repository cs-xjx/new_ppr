% function main_mut()

clear all;
addpath('D:\Polyspace\metis-5.1.0\metismex-dgleich\');
% addpath('D:\Polyspace\metis-5.1.0\metismex-dgleich\yu\');
addpath('D:\Polyspace\metis-5.1.0\metismex-dgleich\amortized\mutual_exclusion\');
addpath('D:\Polyspace\metis-5.1.0\metismex-dgleich\test\');
% addpath('D:\Polyspace\metis-5.1.0\metismex-dgleich\BEAR-2.0\');

% get a
fpath = 'D:\Polyspace\metis-5.1.0\metismex-dgleich\dataset\';
%fpath = 'E:\datasets\Matlab\';
% ds = '494_bus';
% ds = 'ego-Facebook';
% ds = 'wiki-Vote';
% ds = 'p2p-Gnutella31';
% ds = 'email-EuAll';
% ds = 'web-NotreDame';
ds = 'web-BerkStan';
% ds = 'web-Google';
% ds = 'com-Youtube';
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
src = [267061 394151 115798 510286];
tar = {[685279 685246 655173], [685259 562714 685236], [501402 308993 203056], [685313 685275 534636]};
qu = {[454904 644168 267061 510286 310124 283143 564815 115798 76720 394151], ...
    [667497 313129 152236 158152 314017 424609 564501 522972 529530 394151], ...
    [369398 140008 267061 510286 671169 446267 339572 115798 430582 517411], ...
    [208392 461934 267061 510286 179972 172815 273361 428037 600419 394151], ...
    [114825 321962 324594 510286 382904 280498 174924 115798 140691 550351]};

% a = sparse([
% 0, 0, 0, 0, 1, 0, 0, 0, 0;
% 0, 0, 1, 1, 0, 0, 0, 0, 0;
% 0, 0, 0, 0, 0, 1, 0, 0, 0;
% 0, 0, 0, 0, 1, 0, 0, 0, 1;
% 0, 1, 0, 0, 0, 0, 1, 0, 0;
% 0, 1, 0, 0, 0, 0, 0, 0, 0;
% 0, 0, 1, 0, 0, 0, 0, 0, 0;
% 0, 0, 0, 1, 0, 0, 0, 0, 0;
% 1, 0, 0, 0, 0, 0, 0, 0, 0
%     ]);
% src = [1 2 1 2 5];
% tar = {[2 7 14 17], [119 5 6 17], [4 9 6], [90 5], [9 117]};
% qu = {[1 2 6]};
% src = [1 2 5];
% tar = {[4 7 2 4], [3 5 6 9], [9 2]};
% qu = {[1 2 6], [3 7]};

% a = sparse([
% 0, 1, 0, 0, 0, 0, 0;
% 0, 0, 0, 0, 0, 0, 0;
% 0, 0, 0, 0, 0, 0, 0;
% 0, 0, 0, 0, 0, 0, 0;
% 0, 0, 0, 0, 0, 0, 0;
% 0, 0, 0, 0, 0, 0, 0;
% 0, 0, 0, 0, 0, 0, 0
%     ]);
% src = [1 1];
% tar = {[3 4], [5 6 7]};
% qu = {[1 2 6]};

% set c, kmax, num_blks, low rank
c = 0.6;            % decay factor
kmax = 20;         % # of iterations
n_blocks = 500;        % # of partitions
r = 20;        % low rank of svds


%qu_len = 100;       % query length
%qu_num = 10;        % # of queries
%ue_num = 4;  % 8;   % # of uncertain edges (source nodes)
%ue_deg = 10; % 4;   % # of uncertain target nodes for each uncertain edge

uede = ['-', num2str(numel(src)), '_', num2str(numel(tar{1}) + 1)];
% algo
% [res1] = exhPPR_MutualExclusion(a, c, kmax, src, tar, qu, ds);
% save(['exh_mut-', ds, uede], 'res1');

standard_res  = load(['exh_mut-', ds, uede]);
res1 = standard_res.res1;

% [res2] = exhPPR(a, c, kmax, src, tar, qu, ds);

% [res_ea] = exhApxPPR_MutualExclusion(a, c, n_blocks, src, tar, qu, r, ds);
% save(['eap_mut-', num2str(n_blocks), ds, uede], 'res_ea');

% [res_cp] = collPPR_MutualExclusion(a, c, kmax, src, tar, qu, ds);
% save(['cp_mut-', ds, uede], 'res_cp');

% [res_cap] = collApxPPR_MutualExclusion(a, c, n_blocks, src, tar, qu, r, ds);
% save(['cap_mut-', num2str(n_blocks), ds, uede], 'res_cap');

% [res_ca2p] = collApx2PPR_MutualExclusion(a, c, n_blocks, src, tar, qu, r, ds);
% save(['ca2p_mut-', ds, uede], 'res_ca2p');

% [res_fp] = flatPPR_MutualExclusion(a, c, kmax, src, tar, qu, ds);
% save(['fp_mut-', ds, uede], 'res_fp');

[res_fap] = flatApxPPR_MutualExclusion(a, c, n_blocks, src, tar, qu, r, ds);
save(['fap_mut-', num2str(n_blocks), ds, uede], 'res_fap');

% [res_bear] = BearS_uncertain_mut(a, c, src, tar, qu, ds);
% save(['bear_mut-', ds, uede], 'res_bear');

% [res_uppr] = UPPR_MutualExclusion(a, c, src, tar, qu, n_blocks,ds);
% save(['uppr_mut-', num2str(n_blocks), ds, uede], 'res_uppr');

% [res3] = new_ppr_mut_inc(a, c, kmax, src, tar, qu, ds);
% [res3] = new_ppr_mut_nonzerotar(a, c, kmax, src, tar, qu, ds);

% check accuracy
top_k = 50;
[sum_err, mAP, nDCG, kendall, spearman] = compute_precision(res1, res_fap, top_k, 1);
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

% qu_num = numel(qu);
% sum_err = 0;
% for qi = 1 : qu_num
%     sum_err = sum_err + norm(res1{qi} - res_ea{qi});
% end
% fprintf('\n\n========== accuracy ============\n');
% sum_err
% fprintf(' >>         error              :  %f   \n\n\n', sum_err);



