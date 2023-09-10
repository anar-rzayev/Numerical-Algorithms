% Use the power method to approximate the dominant eigenvalue and an associated eigenvector of A given a nonzero vector x(0)
% Terminate when $\frac{\left\|x^{(k)}-x^{(k-1)}\right\|_{\infty}}{\| x^{(k)}||_{\infty}}<\epsilon$ or $k \geq N$, and output $x^{(k)}$


function [mu x iter] = power_method(A, x0, epsilon, N)
    x = x0;
    function [a] = helper(xp)
        desiredVal = norm(xp, Inf);
        a = 1;
        while abs(xp(a)) ~= desiredVal
            a = a + 1;
        end
    end
    iter = 1;
    x = x / norm(x, inf);
    while iter <= N
        y = A * x;
        mu = norm(y, inf);
        err = norm(x - (y / mu), Inf) / norm(x, Inf);
        x = y / mu;
        if err < epsilon
            break
        end
        iter = iter + 1;
    end

end
