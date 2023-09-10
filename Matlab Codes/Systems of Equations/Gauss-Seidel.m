% Use Gauss-Seidel method to solve Ax = b given an initial approximation x(0)


function [xg Ng] = GS(A, b, x0, epsilon, N)
    [n, n] = size(A);
    XO = x0;
    xg = x0;
    Ng = 0;
    for g = 1:N
        x = zeros(n, 1);
        for i = 1 : n      
            x(i) = (1/A(i,i)) * (-A(i, 1:(i-1)) * x(1:(i-1)) - A(i, (i+1):n) * XO((i+1):n) + b(i));
        end
        if norm(x - XO, inf) / norm(x, inf) < epsilon
            xg = x
            Ng = g;
            break
        else
            XO = x;
            xg = XO;
            Ng = Ng + 1;
        end
    end
end