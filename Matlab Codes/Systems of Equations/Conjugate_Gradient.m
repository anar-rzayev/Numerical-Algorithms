% Use the conjugate gradient method to solve Ax = b given an initial approximation x(0)
% Terminate when $\frac{\left\|x^{(k)}-x^{(k-1)}\right\|_{\infty}}{\| x^{(k)}||_{\infty}}<\epsilon$ or $k \geq N$, and output $x^{(k)}$


function [xc Nc] = conjugate_gradient(A, b, x0, epsilon, N)
    r0 = b - A * x0;
    v1 = r0;
    x_tracker = x0;
    r_tracker = r0;
    v_tracker = v1;
    xc = x0;
    Nc = 0;
    for c = 1:N
        r_tracker
        t_k = dot(r_tracker, r_tracker) / dot(v_tracker, A * v_tracker);
        xc = x_tracker + t_k * v_tracker;
        r_tracker_1 = r_tracker - t_k * (A * v_tracker);
        s_k = dot(r_tracker_1, r_tracker_1) / dot(r_tracker, r_tracker);
        v_tracker = r_tracker_1 + s_k * v_tracker;
        r_tracker = r_tracker_1
        if norm(xc - x_tracker, inf) / norm(xc, inf) < epsilon
            Nc = c;
            break
        else
            x_tracker = xc;
            Nc = Nc + 1;
        end
    end

end