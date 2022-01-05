% Given the following pairs of data
% (x0, f(x0)), (x1, f(x1)), ... , (xn, f(xn)), wherew x0 < x1 < ... < xn
% find the numbers a_j, b_j, c_j, d_j of the cubic spline interpolant S.
% For x_j <= x <= x_{j+1}, S(x) = S_j(x) = a_j + b_j * (x - x_j) + c_j * (x - x_j) ** 2 + d_j * (x - x_j) ** 3


function [a b c d] = cubic_spline(x, f, opt, g0, gn)
    if opt == 0
        n = length(f) - 1;
        h = zeros(n, 1);
        alpha = zeros(n - 1, 1);
        for i = 1 : n 
            h(i) = x(i + 1) - x(i);
        end
        
        for j = 1 : (n - 1)
            alpha(j) = (3 / h(j + 1)) * (f(j + 2) - f(j + 1)) - (3 / h(j)) * (f(j + 1) - f(j))
            
        end
        
        l = ones(n + 1, 1);
        mu = zeros(n + 1, 1);
        z = zeros(n + 1, 1);
        
        for k = 1 : (n - 1)
            l(k + 1) = 2 * (x(k + 2) - x(k)) - h(k) * mu(k);
            mu(k + 1) = h(k + 1) / l(k + 1);
            z(k + 1) = (alpha(k) - h(k) * z(k)) / l(k + 1);
        end
        
        c = zeros(n + 1, 1);
        d = zeros(n, 1);
        b = zeros(n, 1);
        a = f(1:n);
        
        for t = n : -1 : 1
            c(t) = z(t) - mu(t) * c(t + 1);
            b(t) = (f(t + 1) - f(t)) / h(t) - h(t) * (c(t + 1) + 2 * c(t)) / 3;
            d(t) = (c(t + 1) - c(t)) / (3 * h(t));
        end
        
        c = c(1:n);
        
    else 
        n = length(f) - 1;
        h = zeros(n, 1);
        alpha = zeros(n + 1, 1);
        for i = 1 : n 
            h(i) = x(i + 1) - x(i);
        end
        
        alpha(1) = 3 * (f(2) - f(1))/ h(1) - 3 * g0;
        alpha(n + 1) = 3 * gn - 3 * (f(n + 1) - f(n)) / h(n);
        
        for j = 1 : (n - 1)
            alpha(j + 1) = (3 / h(j + 1)) * (f(j + 2) - f(j + 1)) - (3 / h(j)) * (f(j + 1) - f(j))
        end
        
        l = ones(n + 1, 1);
        mu = zeros(n + 1, 1);
        z = zeros(n + 1, 1);
        
        l(1) = 2 * h(1);
        mu(1) = 0.5;
        z(1) = alpha(1) / l(1);
        
        for k = 1 : (n - 1)
            l(k + 1) = 2 * (x(k + 2) - x(k)) - h(k) * mu(k);
            mu(k + 1) = h(k + 1) / l(k + 1);
            z(k + 1) = (alpha(k + 1) - h(k) * z(k)) / l(k + 1);
        end
        
        l(n + 1) = h(n) * (2 - mu(n));
        z(n + 1) = (alpha(n + 1) - h(n) * z(n)) / l(n + 1);
        c = zeros(n + 1, 1);
        c(n + 1) = z(n + 1);
        d = zeros(n, 1);
        b = zeros(n, 1);
        a = f(1:n);
        
        for t = n : -1 : 1
            c(t) = z(t) - mu(t) * c(t + 1);
            b(t) = (f(t + 1) - f(t)) / h(t) - h(t) * (c(t + 1) + 2 * c(t)) / 3;
            d(t) = (c(t + 1) - c(t)) / (3 * h(t));
        end
        
        c = c(1:n);
    end
end

