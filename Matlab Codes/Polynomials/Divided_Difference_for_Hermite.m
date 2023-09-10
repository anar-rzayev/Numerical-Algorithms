% Given the following pairs of data
% (x0, f(x0), f`(x0)), (x1, f(x1), f`(x1)), ... , (xn, f(xn), f`(xn))
% find the numbers $a_{k}$ for the $n$th interpolating polynomial $P_{n}(x)=a_{0}+a_{1}\left(x-x_{0}\right)+a_{2}\left(x-x_{0}\right)\left(x-x_{1}\right)+\cdots+a_{n}\left(x-x_{0}\right) \cdots\left(x-x_{n-1}\right)$, 
% and the numbers $b_{k}$ for the $n$th interpolating polynomial $P_{n}(x)=b_{0}+b_{1}\left(x-x_{n}\right)+b_{2}\left(x-x_{n}\right)\left(x-x_{n-1}\right)+\cdots+b_{n}\left(x-x_{n}\right) \cdots\left(x-x_{1}\right)$.


function [a b] = divided_diff_for_Hermite(x, f, g)
    n = length(f) - 1;
    a = zeros(2 * n + 2, 1);
    b = zeros(2 * n + 2, 1);
    matrix = zeros(2 * n + 2, 2 * n + 2);
    z = zeros(2 * n + 2, 1);
    for i = 1 : (n + 1)
        z(2 * i - 1) = x(i);
        z(2 * i) = x(i);
        matrix(2 * i - 1, 1) = f(i);
        matrix(2 * i, 1) = f(i);
        matrix(2 * i, 2) = g(i);
        if i ~= 1
            matrix(2 * i - 1, 2) = (matrix(2 * i - 1, 1) - matrix(2 * i - 2, 1)) / (z(2 * i - 1) - z(2 * i - 2));
        end
    end
    for t = 3 : (2 * n + 2)
        for j = 3 : t
            matrix(t, j) = (matrix(t, j - 1) - matrix(t - 1, j - 1)) / (z(t) - z(t - j + 1));
        end
    end
    for k = 1 : (2 * n + 2)
        a(k) = matrix(k, k);
    end
    b = transpose(matrix(2 * n + 2, :));
end