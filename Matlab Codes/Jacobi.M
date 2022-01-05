% Use Jacobi method to solve Ax = b given an initial approximation x(0)


function [xj Nj] = Jacobi(A, b, x0, epsilon, N)
    D = diag(diag(A));
    T = inv(D) * (D - A);
    new_c = inv(D) * b;
    tracker = x0;
    xj = x0;
    Nj = 0;
    for j = 1:N
        xj = T *i tracker + new_c;
        if norm(xj - tracker, inf) / norm(xj, inf) < epsilon
            Nj = j;
            break
        else
            tracker = xj;
            Nj = Nj + 1;
        end
    end
end