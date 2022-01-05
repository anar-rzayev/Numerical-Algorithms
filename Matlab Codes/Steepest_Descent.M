% Use the steepest descent method to solve Ax = b given an initial approximation x(0)
% Terminate when $\frac{\left\|x^{(k)}-x^{(k-1)}\right\|_{\infty}}{\| x^{(k)}||_{\infty}}<\epsilon$ or $k \geq N$, and output $x^{(k)}$


function [xs Ns] = steepest_descent(A, b, x0, epsilon, N)
    tracker = x0;
    xs = x0;
    Ns = 0;
    for s = 1:N
        t_k = dot(b - A * tracker, b - A * tracker) / dot(b - A * tracker, A * (b - A * tracker));
        xs = tracker + t_k * (b - A * tracker);
        if norm(xs - tracker, inf) / norm(xs, inf) < epsilon
            Ns = s;
            break
        else
            tracker = xs;
            Ns = Ns + 1;
        end
    end
end