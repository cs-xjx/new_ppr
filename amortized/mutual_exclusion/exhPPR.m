function [res] = exhPPR(a, c, kmax, src, tar, qu, ds) 
time_total = tic;
fprintf('\n >> Start exhPPR ... \n');

%fpath = 'D:\MyDoc\Weiren\datasets\matlab\';
% fpath = 'E:\datasets\Matlab\';
% fname = 'soc-LiveJournal1';
% 
% fn = [fpath, fname, '.mat'];

% datasets available: 
%
% https://www.cise.ufl.edu/research/sparse/matrices/SNAP/
%
% as-735 dblp_03-05 dblp_06-08 dblp_09-11
% cit-HepPh p2p-Gnutella08 email-Enron cit-Patents web-Stanford 
% web-BerkStan web-Google email-EuAll uk-2002 soc-LiveJournal1


% c = 0.8;            % decay factor
% kmax = 100;         % # of iterations
% npart = 500;        % # of partitions
% qu_len = 100;       % query length
% qu_num = 10;        % # of queries
% ue_num = 4;  % 8;   % # of uncertain edges (source nodes)
% ue_deg = 10; % 4;   % # of uncertain target nodes for each uncertain edge


%%%%%%%%%%%%%%%%%%%% define a, qu %%%%%%%%%%%%%%%%%%%%
%
% load(fn);
% a = Problem.A;
% [~, idx] = maxk(sum(a), min(size(a,1), 3*qu_len));  % top-(qu_len) query nodes
% qu = cell(qu_num,1);
% for i=1:qu_num
%     qu{i} = idx(randperm(numel(idx), qu_len));
% end

% fpath = 'D:\MyDoc\Weiren\datasets\matlab\';
% %fpath = 'E:\datasets\Matlab\';
% ds = 'as-735';
% fname = [fpath, ds, '.mat'];
% load(fname);
% a = Problem.A;
% n = size(a,1);
% m = nnz(a);
% 
% src = [1912 3059];
% tar = {[3356 2221 1154], [1534 306 4071]};
% qu = {[260 553], [2658 145], [1163 3525], [2930 3248], [1457 1100]};

% a = sparse([
% 0	1	0	0	1
% 0	0	1	1	0
% 0	0	0	0	1
% 0	0	1	0	0
% 0	0	0	0	0
% ]);
% 
% src = [1,4,5];
% %tar = {[2,3,5],[3,6,7],[8,9]};
% tar = {[2,3,5],[3,4],[2,5]};
% 
% %qu = {[4 5],[2 3 5],[1 3]};
% qu = {[3 5]};

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
% tar = {[2 7 14 17], [119 8 6 17], [4 9 6], [90 5], [9 117]};
% qu = {[1 2 6]};
% 
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % % common parameter 
% c = 0.8;            % decay factor
% kmax = 100;         % # of maximum iterations
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% qu_num = numel(qu);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%d = sum(a,2);

%%%%%%%%%%%%%%%%% define src, tar %%%%%%%%%%%%%%%%%%%
%
% [~, src_idx] = maxk(d, min(ue_num*3, na));
% src = sort(src_idx(randperm(numel(src_idx), ue_num)));
% nsrc = numel(src);
% tar = cell(nsrc,1);
% for i=1:nsrc
%     tar_idx = find(a(src(i),:)==0);
%     tar{i} = randperm(numel(tar_idx), ue_deg);
% end



%src = [1];
%tar = {[2,3,5]};

% src = [5];
% tar = {[1]};

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

n_certain = size(a,1);
m_certain = nnz(a);

na = size(a,1);
n = max([na, tar{:}]);
a(n,n) = 0;

nsrc = numel(src);
tar0 = cellfun(@(x) [0 x], tar, 'UniformOutput', false);
tars = cell(1,nsrc);
[tars{end:-1:1}] = ndgrid(tar0{end:-1:1});
tars = cat(nsrc+1, tars{:});
tars = reshape(tars,[],nsrc);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

npw = size(tars,1);
qu_num = numel(qu);
res = cell(qu_num,1);

fprintf('\n== BEGIN online query ==\n');
for ii=1:qu_num
    ppr = sparse(n,1);
    if (mod(ii, 1)==0)
        fprintf(' query #%d\n', ii);
    end
    
    for i = 1:npw
        if (mod(i, 500)==0)
            fprintf(' %% of possible worlds: %2.2f\n', i/npw*100);
        end
        
        tar_i = tars(i,:);
        idx = find(tar_i);
        tar_i = tar_i(idx)';
        src_i = src(idx)';
        nsrc_i = numel(tar_i);
        
        if (nsrc_i>0)
            da = sparse(src_i, tar_i, 1, n, n);
        else
            da = sparse(n, n);
        end
        pa = a + da;
        d = 1./sum(pa,2);
        d(~isfinite(d)) = 0;
        pq = pa'*spdiags(d,0,n,n);
        clear da d;
        ppr = ppr + ppr_iter1(pq, qu{ii}, c, kmax);      
    end
    res{ii} = ppr/npw;
end



fprintf('== END online query ==\n');
period_total = toc(time_total);

fprintf('\n\n========== total stats (exhPPR) ============\n\n');

fprintf(' >>         dataset              :  %s \n', ds);
fprintf(' >>      # of nodes              :  %d \n', n_certain);
fprintf(' >>      # of edges              :  %d \n\n', m_certain);

fprintf(' >>  total CPU time              :  %fs \n', period_total);
fprintf(' >>    # of possible worlds      :  %d  \n', npw);
fprintf(' >>    # of queries              :  %d  \n', qu_num);
fprintf(' >>    ave time per query        :  %fs \n', period_total/qu_num);
fprintf(' >>    ave time per query per pw :  %fs \n', period_total/qu_num/npw);
