% Consider a root-finding problem: sinx - e^(-x) = 0
% Use Newton's method for n>=1, p_n = p_{n-1} - f(p_{n-1})/f`(p_{n-1}) to find an approximate solution 

function sol = newton(p0, N, eps)
    tracker = p0;
    for i = 1:N 
        mediator = tracker - (sin(tracker) - exp(-tracker))/(cos(tracker) + exp(-tracker));
        if abs(mediator - tracker)/abs(mediator) < eps
            sol = [tracker; i];
            break;
        else
            tracker = mediator;
        end
    sol = [tracker; N];      
end

