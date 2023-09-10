% Use the SOR method with the optimal choice of w to solve Ax = b given an initial approximation x(0)
% Terminate when $\frac{\left\|x^{(k)}-x^{(k-1)}\right\|_{\infty}}{\| x^{(k)}||_{\infty}}<\epsilon$ or $k \geq N$, or output $x^{(k)}$


function [ws xs Ns] = SOR(A, b, x0, epsilon, N)
    [n, n] = size(A);
    D = diag(diag(A));
    T = inv(D) * (D - A);
    p_T_j = max(abs(eig(T)));
    ws = 2 / (1 + sqrt(1 - p_T_j * p_T_j));
    XO = x0;
    xs = x0;
    Ns = 0;
    for g = 1:N
        x = zeros(n, 1);
        for i = 1 : n      
            x(i) = (1 - ws) * XO(i) + (ws / A(i,i)) * (-A(i, 1:(i-1)) * x(1:(i-1)) - A(i, (i+1):n) * XO((i+1):n) + b(i));
        end
        if norm(x - XO, inf) / norm(x, inf) < epsilon
            xs = x
            Ns = g;
            break
        else
            XO = x;
            xs = XO;
            Ns = Ns + 1;
        end
    end
end