% Consider the nonlinear system: 3 * x(1) - cos(x(2) * x(3)) - 1/2 = 0; 4 * x(1) * x(1) - 625 * x(2) * x(2) + 2 * x(2) - 1 = 0; exp(-x(1) * x(2)) + 20 * x(3) + (10 * pi - 3) / 3 = 0
% Use Newton's method to find an approximate solution
% Terminate when $\frac{\left\|x^{(k)}-x^{(k-1)}\right\|_{\infty}}{\| x^{(k)}||_{\infty}}<\epsilon$ or $k \geq N$, and output $x^{(k)}$


function sol = newton(x0, N, epsilon)
    function F = orig_func(x)
        F = [3 * x(1) - cos(x(2) * x(3)) - 1/2; 4 * x(1) * x(1) - 625 * x(2) * x(2) + 2 * x(2) - 1; exp(-x(1) * x(2)) + 20 * x(3) + (10 * pi - 3) / 3];
    end
    function J = jacobian(x)
        J = [3, sin(x(2) * x(3)) * x(3), sin(x(2) * x(3)) * x(2); 8 * x(1), -1250 * x(2) + 2, 0; exp(-x(1) * x(2)) * (-x(2)), exp(-x(1) * x(2)) * (-x(1)), 20];
    end
    sol = [x0; 0];
    for i = 1:N
        [k, t] = size(sol);
        x_i = sol(1: k - 1);
        F_i = orig_func(x_i);
        J_i = jacobian(x_i);
        y_i = linsolve(J_i, -F_i);
        new_tracker = x_i + y_i;
        
        sol(k) = i;
        sol(1: k - 1) = new_tracker;
        if norm(y_i, inf) / norm(new_tracker, inf) < epsilon
            break
        end
    end
end