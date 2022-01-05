% Given the following pairs of data
% (x0, f(x0)), (x1, f(x1)), ... , (xn, f(xn))
% find the numbers $a_{k}$ for the $n$th interpolating polynomial $P_{n}(x)=a_{0}+a_{1}\left(x-x_{0}\right)+a_{2}\left(x-x_{0}\right)\left(x-x_{1}\right)+\cdots+a_{n}\left(x-x_{0}\right) \cdots\left(x-x_{n-1}\right)$, 
% and the numbers $b_{k}$ for the $n$th interpolating polynomial $P_{n}(x)=b_{0}+b_{1}\left(x-x_{n}\right)+b_{2}\left(x-x_{n}\right)\left(x-x_{n-1}\right)+\cdots+b_{n}\left(x-x_{n}\right) \cdots\left(x-x_{1}\right)$.


function [a b] = newton_divided_diff(x, f)
    n = length(f) - 1;
    a = zeros(n + 1, 1);
    b = zeros(n + 1, 1);
    matrix = zeros(n + 1, n + 1);
    matrix(:, 1) = f;
    for i = 2 : (n + 1)
        for j = 2 : i
            matrix(i, j) = (matrix(i, j - 1) - matrix(i - 1, j - 1)) / (x(i) - x(i - j + 1));
        end
    end
    for k = 1 : (n + 1)
        a(k) = matrix(k, k);
    end
    b = transpose(matrix(n + 1, :));   
end