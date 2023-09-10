% Approximate integral of f(x) from a to b using 8 different Newton-Cotes methods:
% Closed Newton-Cotes n = 1, 2, 3, 4
% Open Newton-Cotes n = 0, 1, 2, 3


function approx = newton_cotes(f, a, b, n, Closed)
    if (Closed == 1)
        if (n == 1)
            x0 = a;
            x1 = b;
            h = b - a;
            approx = (h / 2) * f(x0) + (h / 2) * f(x1);
            
        elseif (n == 2)
            x0 = a;
            x1 = a + (b - a) / 2;
            x2 = b;
            h = (b - a) / 3;
            approx = (h / 3) * f(x0) + (4 * h / 3) * f(x1) + (h / 3) * f(x2);
            
        elseif (n == 3)
            x0 = a;
            x1 = a + (b - a) / 3;
            x2 = a + 2 * (b - a) / 3;
            x3 = b;
            h = (b - a) / 3;
            approx = (3 * h / 8) * f(x0) + (9 * h / 8) * f(x1) + (9 * h / 8) * f(x2) + (3 * h / 8) * f(x3);
            
        else
            x0 = a;
            x1 = a + (b - a) / 4;
            x2 = a + 2 * (b - a) / 4;
            x3 = a + 3 * (b - a) / 4;
            x4 = b;
            h = (b - a) / 4;
            approx = (14 * h / 45) * f(x0) + (64 * h / 45) * f(x1) + (24 * h / 45) * f(x2) + (64 * h / 45) * f(x3) + (14 * h / 45) * f(x4);
            
        end
  
        
    else
        if (n == 0)
            h = (b - a) / 2;
            x0 = a + h;
            approx = 2 * h * f(x0);
            
        elseif (n == 1)
            h = (b - a) / 3;
            x0 = a + h;
            x1 = x0 + h;
            approx = (3 * h / 2) * f(x0) + (3 * h / 2) * f(x1);
            
        elseif (n == 2)
            h = (b - a) / 4;
            x0 = a + h;
            x1 = x0 + h;
            x2 = x1 + h;
            approx = (8 * h / 3) * f(x0) - (4 * h / 3) * f(x1) + (8 * h / 3) * f(x2);
            
        else
            h = (b - a) / 5;
            x0 = a + h;
            x1 = x0 + h;
            x2 = x1 + h;
            x3 = x2 + h;
            approx = (55 * h / 24) * f(x0) + (5 * h / 24) * f(x1) + (5 * h / 24) * f(x2) + (55 * h / 24) * f(x3);
           
        end
        
    end 
    
end

