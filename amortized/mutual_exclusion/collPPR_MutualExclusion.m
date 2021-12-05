function [res] = collPPR_MutualExclusion(a, c, kmax, src, tar, qu, ds) 
time_total = tic;
fprintf('\n >> Start collPPR_MutualExclusion ... \n');

n_certain = size(a,1);
m_certain = nnz(a);

qu_num = numel(qu);
nsrc = numel(src);

n = max([n_certain, tar{:}]);

a(n, n) = 0;

%
% generate possible world (mutual exclusive semantics)
%
time_gen_pwd = tic;
fprintf(' >> Generate Possible World ... \n');
tar0 = cellfun(@(x) [0 x], tar, 'UniformOutput', false);
tars = cell(1,nsrc);
[tars{end:-1:1}] = ndgrid(tar0{end:-1:1});
tars = cat(nsrc+1, tars{:});
tars = reshape(tars,[],nsrc);  % tars: possible world
npw = size(tars, 1);
period_gen_pwd = toc(time_gen_pwd);
fprintf('  # of possible worlds            :  %d  \n', npw);
fprintf('  Time (generate possible worlds) :  %fs \n', period_gen_pwd);
%
%

%
% 1. collapse all pws
%
time_collapse = tic;
fprintf(' >> 1. Collapse All Possible Worlds ... \n');
% qavg: avg of all transition matrices
% qavg = sparse([], [], [], n, n, floor(n * n * 0.1)); % 10% sparse
qavg = sparse(n, n);
% enumerate all possible worlds by loop
for i = 1 : npw
    % find the uncertain part of the pw
    tar_i = tars(i,:);
    idx = find(tar_i);
    tar_i = tar_i(idx);
    src_i = src(idx);
    
    % pq: transition matrix
    pa = a + sparse(src_i, tar_i, 1, n, n);
    d = 1./sum(pa,2);
    d(~isfinite(d)) = 0;
    pq = pa' * spdiags(d,0,n,n);
    clear d pa;
    
    % get the avg of transition matrixes of all pws
    qavg = qavg + pq / npw;
end
period_collapse = toc(time_collapse);

%
% 2. online query 
%
time_online_query = tic;
fprintf(' >> 2. Online Query ... \n');
res = cell(qu_num, 1);
for qi = 1 : qu_num
    fprintf(' query #%d\n', qi);
    res{qi} = ppr_iter(c, qavg, qu{qi}, kmax);
end
period_online_query = toc(time_online_query);

period_total = toc(time_total);

fprintf('\n========== collPPR_MutualExclusion ============\n');

fprintf(' >> Dataset                                 :  %s \n',   ds);
fprintf(' >      # of nodes (certain)                :  %d \n',   n_certain);
fprintf(' >      # of edges (certain)                :  %d \n',   m_certain);
fprintf(' >      # of possible worlds                :  %d \n',   npw);
fprintf(' >      # of queries                        :  %d \n\n', qu_num);

fprintf(' >> Parameters\n');
fprintf(' >      decay factor (c)                    :  %d \n',   c);
fprintf(' >      # of iterations                     :  %d \n', kmax);
% fprintf(' >      # of partitions                     :  %d \n', n_blocks);
% fprintf(' >      low-rank SVD (r)                    :  %d \n',   r);
fprintf('\n');

fprintf(' >> Total Time                              :  %fs \n', period_total);
fprintf(' >    1) Preprocessing - Collapse           :  %f s \n', period_collapse);
fprintf(' >    2) Online Query (# of queries = %d)   :  %f s \n\n', qu_num, period_online_query);

fprintf(' >> Average Amortised Total Time per Query              :  %fs \n', period_total / qu_num);
fprintf(' >> Average Amortised Total Time Time per Query per PW  :  %fs \n', period_total / qu_num / npw);
fprintf(' >> Average Amortised Preprocessing Time per Query      :  %f s \n', period_collapse / qu_num);


whos