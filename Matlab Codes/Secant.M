% Consider a root-finding problem: sinx - e^(-x) = 0
% Use the secant method for n>=2, p_n = p_{n-1} - f(p_{n-1}) * (p_{n-1} - p_{n-2})/f(p_{n-1}) - f(p_{n-2}) to find an approximate solution 


function sol = secant(p0, p1, N, eps)
    tracker_1 = p0;
    tracker_2 = p1;
    for i = 1:N 
        fun_tracker_1 = sin(tracker_1) - exp(-tracker_1);
        fun_tracker_2 = sin(tracker_2) - exp(-tracker_2);
        mediator = tracker_2 - fun_tracker_2 * (tracker_2 - tracker_1) / (fun_tracker_2 - fun_tracker_1);
        if abs(mediator - tracker_2)/abs(mediator) < eps
            sol = [mediator; i + 1];
            break;
        else
            tracker_1 = tracker_2;
            tracker_2 = mediator;
        end
    sol = [tracker_2; N + 1];     
end