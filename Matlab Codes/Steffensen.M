% Consider the following fixed-point problem: x = g(x) = sqrt(e^x / 3)
% Use Steffen's method to generate sequence {p_n} to approximate a solution


function sol = steffensen(p0, N, eps)
    tracker = p0;
    for i = 1:N 
        p1 = sqrt(exp(tracker) / 3);
        p2 = sqrt(exp(p1) / 3);
        mediator = tracker - ((p1 - tracker) ^ 2) / (p2 - 2 * p1 + tracker);
        if abs(mediator - p2)/abs(mediator) < eps
            sol = [mediator; i];
            break;
        else
            tracker = mediator
        end
    sol = [tracker; N];     
end