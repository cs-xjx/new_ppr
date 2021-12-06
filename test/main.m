function main()

clear all;
% addpath('D:\Polyspace\metis-5.1.0\metismex-dgleich\');
% addpath('D:\Polyspace\metis-5.1.0\metismex-dgleich\yu\');
addpath('D:\Polyspace\metis-5.1.0\metismex-dgleich\amortized\mutual_exclusion\');
addpath('D:\Polyspace\metis-5.1.0\metismex-dgleich\amortized\');
% addpath('D:\Polyspace\metis-5.1.0\metismex-dgleich\BEAR-2.0\');

% get a
fpath = 'D:\Polyspace\metis-5.1.0\metismex-dgleich\dataset\';
%fpath = 'E:\datasets\Matlab\';
% ds = '494_bus';
ds = 'ego-Facebook';
% ds = 'wiki-Vote';
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
src = [764 3754 2064 1417];
tar = {[740 2531 327], [1962 1800 2050], [3235 1484 2100], [3546 2538 826]};
qu = {[260 553], [2658 145], [1163 3525], [2930 3248], [1457 1100]};

% src = [6105 59906 60051 46510];
% tar = {[34232 60389 62703], [50086 62637 49580], [2235 58449 47419], [41024 62671 62620]};
% qu = {[869 1787], [29957 27856], [38140 49031], [50589 57839], [23886 50383]};

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
% % tar = {[2 7 14 17], [119 5 6 17], [4 9 6], [90 5], [9 117]};
% tar = {[4 7 2], [3 5 6 9], [6 9 4], [1 5], [9 2]};
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
c = 0.8;            % decay factor
kmax = 30;         % # of iterations
n_blocks = 3;        % # of partitions
r = 20;        % low rank of svds


%qu_len = 100;       % query length
%qu_num = 10;        % # of queries
%ue_num = 4;  % 8;   % # of uncertain edges (source nodes)
%ue_deg = 10; % 4;   % # of uncertain target nodes for each uncertain edge

% algo
[res1] = exhPPR_MutualExclusion(a, c, kmax, src, tar, qu, ds);
% [res2] = exhPPR(a, c, kmax, src, tar, qu, ds);
% [res2] = exhApxPPR_MutualExclusion(a, c, n_blocks, src, tar, qu, r, ds);
% [res3] = exhApxPPR2_MutualExclusion(a, c, n_blocks, src, tar, qu, r, ds);
% [res3] = collPPR_MutualExclusion(a, c, kmax, src, tar, qu, ds);
% [res3] = collApxPPR_MutualExclusion(a, c, n_blocks, src, tar, qu, r, ds);
% [res3] = collApx2PPR_MutualExclusion(a, c, n_blocks, src, tar, qu, r, ds);
% [res3] = flatPPR_MutualExclusion(a, c, kmax, src, tar, qu, ds);
% [res3] = flatApxPPR_MutualExclusion(a, c, n_blocks, src, tar, qu, r, ds);
% [res3] = BearS_uncertain_mut(a, c, src, tar, qu, ds);
% [res2] = UPPR_MutualExclusion(a, c, src, tar, qu, n_blocks,ds);
% [res2] = new_ppr_mut(a, c, kmax, src, tar, qu, ds);
% [res3] = new_ppr_mut_nonzerotar(a, c, kmax, src, tar, qu, ds);
[res3] = new_ppr_mut_inc(a, c, kmax, src, tar, qu, ds);

% check accuracy
qu_num = numel(qu);
sum_err = 0;
for qi = 1:qu_num
    sum_err = sum_err + norm(res1{qi} - res3{qi});
end
fprintf('\n\n========== accuracy ============\n');
sum_err
fprintf(' >>         error              :  %f  \n\n\n', sum_err);



