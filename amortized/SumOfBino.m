function U = SumOfBino(X_B, u1)
%SumOfBino  summation of binomial distributions,
%           the principle is in the paper "The distribution of a sum of binomial random variables"
% u1  # of the X_B(:, 1)

n = size(X_B, 1);

n1 = X_B(1, 1);
p1 = X_B(1, 2);
U = zeros(1, u1 + 1);

U(1) = (1 - p1) ^ n1;
for j = 2 : u1 + 1
    if j - 1 <= n1
        U(j) = ((n1 - (j - 1) + 1) / (j - 1)) * (p1 / (1 - p1)) * U(j - 1);
    end
end
if u1 == 1
    return;
end

Z = zeros(1, u1 + 1);
for i = 2 : n
    Y = U;
    
    n2 = X_B(i, 1);
    p2 = X_B(i, 2);
    Z(1) = (1 - p2) ^ n2;
    for  j = 2 : u1 + 1
        if j - 1 <= n2
            Z(j) = ((n2 - (j - 1) + 1) / (j - 1)) * (p2 / (1 - p2)) * Z(j - 1);
        end
    end
    
    for j = 0 : u1
        U(j + 1) = 0;
        for k = 0 : j
            U(j + 1) = U(j + 1) + Y(k + 1) * Z(j - k + 1);
        end
    end
    Z(:) = 0;
end


end

