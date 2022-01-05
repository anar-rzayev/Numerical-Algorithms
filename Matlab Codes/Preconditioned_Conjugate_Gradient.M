% Use the preconditioned conjugate gradient method to solve Ax = b given an initial approximation x(0)
% Terminate when $\frac{\left\|x^{(k)}-x^{(k-1)}\right\|_{\infty}}{\| x^{(k)}||_{\infty}}<\epsilon$ or $k \geq N$, and output $x^{(k)}$
% Use C = D^(1/2) for the preconditioned conjugate gradient method where D is a diagonal matrix and has the diagonal entries the same as A


function [xp Np] = preconditioned_conjugate_gradient(A, b, x0, epsilon, N)
    C = diag(sqrt(diag(A)));
    M = inv(C);
    r0 = b - A * x0;
    w0 = M * r0;
    v1 = M * w0;
    alpha = dot(w0, w0);
    r_tracker = r0;
    w_tracker = w0;
    v_tracker = v1;
    xp = x0;
    Np = 1;
    for p = 1:N
        Np = p + 1
        if norm(v_tracker, inf) < epsilon
            break
        else
            u = A * v_tracker;
            t = alpha / dot(v_tracker, u);
            xp = xp + t * v_tracker;
            r_tracker = r_tracker - t * u;
            w_tracker = M * r_tracker;
            betha = dot(w_tracker, w_tracker);
            if abs(betha) < epsilon & norm(r_tracker, inf) < epsilon
                break
            else
                s = betha / alpha;
                v_tracker = M * w_tracker + s * v_tracker;
                alpha = betha;
            end
        end
    end
end