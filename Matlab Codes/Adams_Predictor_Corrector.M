% Use the Adams 4th-order predictor-corrector method to approximate the solution of the system
% y' = f(t, y), a <= t <= b, y(a) = alpha
% at (N + 1) equally spaced numbers in the interval [a, b]
% Use the Runge-Kutta method of order 4 to calculate the starting values w0, w1, w2, w3

function w = IVP_Adams_Predictor_Corrector(f, a, b, alpha, N)
    M = 4;
    h = (b - a) / N;
    w = ones(N, 1);
    w(1) = alpha;
    
    coeffic = ones(1, M);
    timer = ones(1, M);
    coeffic(1) = alpha;
    timer(1) = a;
    
    
    for i = 1 : 3
        
        K1 = h * f(timer(i), coeffic(i)); 
        K2 = h * f(timer(i) + h / 2, coeffic(i) + K1 / 2);
        K3 = h * f(timer(i) + h / 2, coeffic(i) + K2 / 2);
        K4 = h * f(timer(i) + h, coeffic(i) + K3);
        coeffic(i + 1) = coeffic(i) + (K1 + 2 * K2 + 2 * K3 + K4) / 6;
        timer(i + 1) = a + i * h;
        w(i + 1) = coeffic(i + 1);
        
    end
    
    for j = 4 : N
        
        prediction = coeffic(4) + (h / 24) * (55 * f(timer(4), coeffic(4)) - 59 * f(timer(3), coeffic(3)) + 37 * f(timer(2), coeffic(2)) - 9 * f(timer(1), coeffic(1)));
        real = coeffic(4) + (h / 24) * (9 * f(a + j * h, prediction) + 19 * f(timer(4), coeffic(4)) - 5 * f(timer(3), coeffic(3)) + f(timer(2), coeffic(2)));
        
        for k = 1 : 3
            timer(k) = timer(k + 1);
            coeffic(k) = coeffic(k + 1);
        end
        
        timer(4) = a + j * h;
        coeffic(4) = real;
        w(j + 1) = real;
        
    end
end

