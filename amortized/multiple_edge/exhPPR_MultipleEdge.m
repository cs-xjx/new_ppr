function [res] = exhPPR_MultipleEdge(a, c, kmax, src, tar, qu, ds) 
time_total = tic;
fprintf('\n >> Start exhPPR_MultipleEdge ... \n\n');

n_certain = size(a,1); % n_certain: # of nodes in graph
m_certain = nnz(a);

qu_num = numel(qu);
nsrc = numel(src);

n = max([n_certain, tar{:}]);

a(n, n) = 0;



%
% generate possible world (multiple edge semantics)
%
time_gen_pwd = tic;
fprintf(' >> Generate Possible World ... \n');
tars = cartprod_MultipleEdge(tar, nsrc);
[npw, ~] = size(tars); % numOfWorlds: # of all possible worlds
period_gen_pwd = toc(time_gen_pwd);
fprintf('  # of possible worlds            :  %d  \n', npw);
fprintf('  Time (generate possible worlds) :  %fs \n', period_gen_pwd);
%
%

%
% 1. online query 
%
fprintf(' >> 1. Online Query ... \n');
res = cell(qu_num, 1);
period_online_query = 0;
for qi = 1 : qu_num
    % pSum: sum of all r
    psum = sparse(n, 1);
    fprintf(' query #%d\n', qi);
    % enumerate all possible worlds by loop
    for i = 1 : npw
        pa = a;
        for j = 1 : nsrc
            pa = pa + sparse(src(j), tars{i, j}, 1, n, n);
        end
        % pq: transition matrix
        d = 1./sum(pa,2);
        d(~isfinite(d)) = 0;
        pq = pa' * spdiags(d,0,n,n);
        % get PPR
        time_online_query = tic;
        psum = psum + ppr_iter(c, pq, qu{qi}, kmax);
        period_online_query = period_online_query + toc(time_online_query);
    end
    res{qi} = psum / npw;
end
% period_online_query = toc(time_online_query);
period_total = toc(time_total);

fprintf('\n========== exhPPR_MultipleEdge ============\n');

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
fprintf(' >    1) Online Query (# of queries = %d)   :  %fs \n\n', qu_num, period_online_query);

fprintf(' >> Average Amortised Total Time per Query              :  %fs \n', period_total / qu_num);
fprintf(' >> Average Amortised Total Time Time per Query per PW  :  %fs \n', period_total / qu_num / npw);

whos

