% Compare the various methods of calculation
% a = sparse([1 2 3 5 5 7 7], [2 4 4 4 6 5 6], 1, 8, 8); % graph 1
% src = [3 3 7];
% tar = {[2], [5 6 8], [4 8]};

a = sparse([
0, 0, 0, 0, 1, 0, 0, 0, 0;
0, 0, 1, 1, 0, 0, 0, 0, 0;
0, 0, 0, 0, 0, 1, 0, 0, 0;
0, 0, 0, 0, 1, 0, 0, 0, 1;
0, 1, 0, 0, 0, 0, 1, 0, 0;
0, 1, 0, 0, 0, 0, 0, 0, 0;
0, 0, 1, 0, 0, 0, 0, 0, 0;
0, 0, 0, 1, 0, 0, 0, 0, 0;
1, 0, 0, 0, 0, 0, 0, 0, 0
    ]);
src = [1 1 2 5];
tar = {[4 7 2 8], [6 9 4], [1 5], [9 2]};
qu = [1 2 6];

n = size(a, 1);
c = 0.8;

suma = sum(a, 2);
d = 1 ./ suma;
d(~isfinite(d)) = 0;
Q = a' * spdiags(d, 0, n, n);
En = eye(n);
R = inv(En - c * Q);
s = sparse(qu, 1, 1 / numel(qu), n, 1);
p0 = R * s;

unique_src = unique(src);
u_nsrc = numel(unique_src);
D = diag(suma);

% I1 = [3 7];
% J1 = [5 4];
% I1 = [1 2];
J1 = [4 6 0 9];
L1_low = diag([2 1]);
L1_high = diag([2 0 1]);

%
% 1. pw1, low dimension
%
idx = find(J1);
I1 = src(idx);
u_I1 = unique(I1);
tJ = J1(idx);
D1_low = D(u_I1, u_I1);
% L1_low = diag([2 1]);
M1_low = sparse(tJ, I1, 1, n, n);
M1_low = M1_low(:, u_I1);
% output
RIJ = R(u_I1, :) * M1_low;
H1_low = inv(D1_low - c * RIJ + R(u_I1, u_I1) * L1_low);
h1_low = H1_low * p0(u_I1)
p1_low = p0 + (c * R * M1_low - R(:, u_I1) * L1_low + En(:, u_I1) * L1_low) * h1_low

%
% 2. pw1, high dimension
%
D1_high = D(unique_src, unique_src);
D1_high(2, 2) = 1;
% L1_high = diag([2 0 1]);
M1_high = sparse(J1 + 1, src, 1, n + 1, n);
M1_high = M1_high(:, unique_src);
M1_high(1, :) = [];
% output
RIJ = R(unique_src, :) * M1_high;
H1_high = inv(D1_high - c * RIJ + R(unique_src, unique_src) * L1_high);
h1_high = H1_high * p0(unique_src)
p1_high = p0 + (c * R * M1_high - R(:, unique_src) * L1_high + En(:, unique_src) * L1_high) * h1_high

% % % % % % % % % % % % % % 
J2 = [4 6 1 9];
k = 3;
u_k = find(unique_src == src(k));
dk = D(unique_src(u_k), unique_src(u_k));
L2_low = diag([2 1 1]);
L2 = diag([2 1 1]);
M2 = sparse(J2 + 1, src, 1, n + 1, n);


%
% 3. pw2, compute by (1)
%
idx2 = find(J2);
I2 = src(idx2);
u_I2 = unique(I2);
post_nsrc = numel(u_I2);
M2_low = M2(:, u_I2);
M2_low(1, :) = [];
u = (dk - 1) * sparse(u_k, 1, 1, post_nsrc, 1) - c * R(u_I2, J2(k)) + R(u_I2, src(k));
supplement = (c * R(src(k), :) * M1_low - R(src(k), u_I1) * L1_low) * H1_low;
new_H1 = zeros(post_nsrc, post_nsrc);
t_index = 1 : post_nsrc;
t_index = t_index(t_index ~= u_k);
new_H1(t_index, t_index) = H1_low;
new_H1(u_k, t_index) = supplement;
new_H1(:, u_k) = sparse(u_k, 1, 1, post_nsrc, 1);

y = new_H1 * u;
H2_3 = new_H1 - (1 / (1 + y(u_k))) * y * new_H1(u_k, :);
h2_3 = H2_3 * p0(u_I2)
p2_3 = p0 + (c * R * M2_low - R(:, u_I2) * L2_low + En(:, u_I2) * L2_low) * h2_3




M2 = M2(:, unique_src);
M2(1, :) = [];

%
% 4. pw2, compute by (2)
%
u = (dk - D1_high(u_k, u_k)) * sparse(u_k, 1, 1, u_nsrc, 1) - c * R(unique_src, J2(k)) + R(unique_src, src(k));
y = H1_high * u;
H2_4 = H1_high - (1 / (1 + y(u_k))) * y * H1_high(u_k, :);
h2_4 = H2_4 * p0(unique_src)
p2_4 = p0 + (c * R * M2 - R(:, unique_src) * L2 + En(:, unique_src) * L2) * h2_4



%
% 5. pw2, compute directly
%
D2 = D(unique_src, unique_src);
% output
RIJ = R(unique_src, :) * M2;
H2_5 = inv(D2 - c * RIJ + R(unique_src, unique_src) * L2);
h2_5 = H2_5 * p0(unique_src)
p2_5 = p0 + (c * R * M2 - R(:, unique_src) * L2 + En(:, unique_src) * L2) * h2_5



