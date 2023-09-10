% Approximate integral of f(x) from a to b using 3 different composite numerical integration methods:
% Composite Trapezoidal rule for n subintervals (output: s1 ) 
% Composite Simpson's rule for n subintervals (output: s2) 
% Composite Midpoint rule for n+2 subintervals (output: s3 )


function [s1 s2 s3] = composite_integration(f, a, b, n)
    
    function [simpson] = simpson(f, a, b, n)
        h = (b - a) / n;
        XI0 = f(a) + f(b);
        XI1 = 0;
        XI2 = 0;
        for i  = 1 : (n - 1)
            X = a + i * h;
            if bitget(i, 1)
                XI1 = XI1 + f(X);
            else
                XI2 = XI2 + f(X);
            end
        end
        XI = h * (XI0 + 2 * XI2 + 4 * XI1) / 3;
        simpson = XI;
    end
    
    function [trap] = trapezoidal(f, a, b, n)
        h = (b - a) / n;
        XI0 = f(a) + f(b);
        XI1 = 0;
        for i  = 1 : (n - 1)
            X = a + i * h;
            XI1 = XI1 + f(X);
        end
        XI = (h / 2) * (XI0 + 2 * XI1);
        trap = XI;
    end
    
    function [compMid] = composite(f, a, b, n)
        h = (b - a) / (n + 2);
        XI1 = 0;
        for i  = 0 : (n / 2)
            X = a + (2 * i + 1) * h;
            XI1 = XI1 + f(X);
        end
        XI = 2 * h * XI1;
        compMid = XI;
    end
    
    s1 = trapezoidal(f, a, b, n);
    s2 = simpson(f, a, b, n);
    s3 = composite(f, a, b, n);

    
end

