function [res] = new_ppr_mut(a, c, kmax, src, tar, qu, ds)
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
% tar = {[4 7 2 4], [3 5 6 9], [6 9 4], [1 5], [9 2]};
% qu = [1 2 6];

% time_total = tic;
fprintf('\n >> Start new_ppr_mut ... \n\n');

n_certain = size(a,1); % n_certain: # of nodes in graph
% m_certain = nnz(a);

qu_num = numel(qu);
nsrc = numel(src);

n = max([n_certain, tar{:}]);

% if n > n_certain
    a(n, n) = 0;
% else
%     n = n_certain;
% end

%
% generate possible world (mutual exclusive semantics)
%
time_gen_pwd = tic;
fprintf(' >> Generate Possible World ... \n');
tar0 = cellfun(@(x) [0 x], tar, 'UniformOutput', false);
% tar0 = tar;
tars = cell(1,nsrc);
[tars{end:-1:1}] = ndgrid(tar0{end:-1:1});
tars = cat(nsrc + 1, tars{:});
tars = reshape(tars, [], nsrc);  % tars: possible world
npw = size(tars, 1);
period_gen_pwd = toc(time_gen_pwd);
fprintf('  # of possible worlds            :  %d  \n', npw);
fprintf('  Time (generate possible worlds) :  %fs \n', period_gen_pwd);
%
%

suma = sum(a, 2);
d = 1 ./ suma;
d(~isfinite(d)) = 0;
Q = a' * spdiags(d, 0, n, n); % transition matrix of a
% D = spdiags(suma(src), 0, nsrc, nsrc);    % outdegree of source nodes
D = suma(src);

I = speye(n);
R = inv(I - c * Q);
% RI_ = I(src, :) * R;
% RI_0 = I(src, :);
% RI_ = RI_0;
% for iter = 1 : kmax
%     RI_ = c * RI_ * Q + RI_0;
% end

res = cell(qu_num, 1);
for qi = 1 : qu_num
    s = sparse(qu{qi}, 1, 1 / numel(qu{qi}), n, 1);
    p0 = R * s;
    
sum_h = zeros(n, 1);
sum_h2 = sparse(n, 1);
% sum_h = sparse(n, n);
% sum_h2 = sparse(n, n);
for i = 1 : npw
    tar_i = tars(i,:);
    idx = find(tar_i);
    tar_i = tar_i(idx);
    src_i = src(idx);
%     D_i = diag(D(src_i));
    D_i = diag(D(idx));
    
%     h = (D_i - c * R(src_i, tar_i) + R(src_i, src_i)) \ I(src_i, :);
    h = (D_i - c * R(src_i, tar_i) + R(src_i, src_i)) \ p0(src_i);
    sum_h = sum_h + I(:, src_i) * h;
%     sum_h = sum_h + h;
    sum_h2 = sum_h2 + I(:, tar_i) * h;
end
sum_x0 = c * sum_h2 - sum_h;
% sum_x0 = c * sum_h2 - I(:, src_i) * sum_h;

sum_xk = iter_compute(c, Q, sum_x0, kmax);
% res{qi} = (p0 + sum_xk + I(:, src_i) * sum_h) / npw;
res{qi} = p0 + sum_xk / npw + sum_h / npw;
res{qi} = (1 - c) * res{qi};
res{qi};

end
% xk = R * x0;
% sum_xk = iter_compute(c, Q, sum_x0, kmax);
% 
% temp_p = (I + sum_xk + sum_h) / npw;
% 
% res = cell(qu_num, 1);
% for qi = 1 : qu_num
% %     s = sparse(qu{qi}, 1, 1 / numel(qu{qi}), n, 1);
% %     p0 = R * s;
%     p0 = sum(R(:, qu{qi}), 2) / numel(qu{qi});
%     res{qi} = temp_p * p0;
%     res{qi}
% end


% period_total = toc(time_total);





