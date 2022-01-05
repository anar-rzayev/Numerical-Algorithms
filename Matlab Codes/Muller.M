% Consider a root-finding problem: f(x) = x^4 + 2.4 * x^3 - 12.95 * x^2 - 34.608 * x + 91.296 = 0
% Use Muller's method for n>=1, to find an approximate solution 


function sol = muller(p0, p1, p2, N, eps)
    function eval = func(x)
        eval = x^4 + 2.4 * x^3 - 12.95 * x^2 - 34.608 * x + 91.296;
    end
    sol = [p0; 1];
    for i = 3:N 
        h_1 = p1 - p0;
        h_2 = p2 - p1;
        sigma_1 = (func(p1) - func(p0)) / h_1;
        sigma_2 = (func(p2) - func(p1)) / h_2;
        d = (sigma_2 - sigma_1) / (h_2 + h_1);
        b = sigma_2 + h_2 * d;
        D = sqrt(b * b - 4 * func(p2) * d);
        if abs(b - D) < abs(b + D)
            E = b + D;
        else
            E = b - D;
        end
        h = -2 * func(p2) / E;
        p = p2 + h;
        sol = [p; i];
        if abs(p - p2) / abs(p) < eps
            break;
        else
            p0 = p1; 
            p1 = p2;
            p2 = p;
        end
        
    end
    
end