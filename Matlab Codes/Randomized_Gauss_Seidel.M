% The version of Gauss-Seidel method studied in class uses a cyclic order. In other words, at each iteration $k$, each element $x_{i}$ is sequentially updated from $i=1$ to $i=n$. Develop a 
% randomized Gauss-Seidel, which shuffles the order of $\{1,2, \ldots, n\}$ at each iteration $k$.
% In specific, at each $k$, sample a permutation of $\{1,2, \ldots, n\}$ uniformly at random using the function randperm, and use such order.
% (Don't use randperm more than required, so that the randperm generates the permutation the same as the reference code. Note that randperm is not truely random, and rng(1) is set before each run of the method, so that you receive same random(?) numbers.)
% Terminate when $\frac{\left\|x^{(k)}-x^{(k-1)}\right\|_{\infty}}{\| x^{(k)}||_{\infty}}<\epsilon$ or $k \geq N$, and output $x^{(k)}$


function [xg Ng] = Randomized_GS(A, b, x0, epsilon, N)
    [n, n] = size(A);
    XO = x0;
    xg = x0;
    Ng = 0;
    for g = 1:N
        x = zeros(n, 1);
        r = randperm(n);
        for i = 1 : n      
            x(r(i)) = (1/A(r(i),r(i))) * (-A(r(i), r(1:(i-1))) * x(r(1:(i-1))) - A(r(i), r((i+1):n)) * XO(r((i+1):n)) + b(r(i)));
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