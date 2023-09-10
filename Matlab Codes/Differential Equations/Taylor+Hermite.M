% Use Taylor's method of order 2 to approximate the solution of the initial-value problem
% y' = f(t, y), a <= t <= b, y(a) = alpha
% at (N + 1) equally spaced numbers in the interval [a, b]
% Then, use cubic Hermite polynomial interpolation to approximate y(t0) for some t0 in [a, b]

function w = IVP_Taylor_Hermite(f, fp, a, b, alpha, N, t0)
    h = (b-a) / N;
    whelp = alpha;
    t = a;
    for k = 1 : N
        f_1 = f(t, whelp);
        fp_1 = fp(t, whelp);
        T = f_1 + h/2 * fp_1;
        inc = t + h;
        if (inc > t0) && (t0 >= t)
            whelper = whelp + h*T;
            thelper = inc;
            f_2 = f(thelper, whelper);
            fp_2 = fp(thelper, whelper);
            break
        end
        whelp = whelp + h * T;
        t = inc;
    end
    t_0 = t0 - t;
    t_1 = t0 - thelper;
    w_1 = whelper - whelp;
    w = whelp + t_0 * f_1 + t_0 * t_0 * ((w_1 / h - f_1) / h) + t_0 * t_0 * t_1 * (((f_2 - w_1 / h) / h - (w_1 / h - f_1) / h) / h);
end

