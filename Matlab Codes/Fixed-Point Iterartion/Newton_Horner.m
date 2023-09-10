% Consider the following polynomial: P(x) = x^4 - 4 * x^2 - 3 * x + 5
% Use Newton's method with Horner's method to find an approximate solution 


function sol = newton_horner(p0, N, eps)
    function [f_q, y, z] = Horner(f, n, x_0)
        f_q = zeros(n, 1);
        f_q(1) = f(1);
        z = f(1);
        for j = 2:n
            f_q(j) = x_0 * f_q(j - 1) + f(j);
            z = x_0 * z + f_q(j);
        end
        y = x_0 * f_q(n) + f(n + 1);
    end
    tracker = p0;
    sol = zeros(6, 1);
    for i = 1:N 
        [Q, p_x_0, p_x_0_prime] = Horner([1;0;-4;-3;5], 4, tracker);
        mediator = tracker - p_x_0 / p_x_0_prime;
        sol = [tracker; i; Q(1); Q(2); Q(3); Q(4)];
        if abs(mediator - tracker)/abs(mediator) < eps
            break;
        else
            tracker = mediator;
        end
    end
end