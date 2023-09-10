N = 600;

X_1 = zeros(1, N - 99);
Y_1 = zeros(1, N - 99);

X_2 = zeros(1, N - 99);
Y_2 = zeros(1, N - 99);

X_3 = zeros(1, N - 99);
Y_3 = zeros(1, N - 99);

for m=100:N
    X_1(m-99) = m; X_2(m-99) = m; X_3(m-99) = m;
    
    A = randn(m,m);
    norm_A = norm(A, "fro");

    [QM,RM]=mgs(A);
    [QH,RH]=qr(A);
    [U, D, V] = svd(A);
    
    Y_1(m-99) = norm(A - QM*RM, "fro") / norm_A;
    Y_2(m-99) = norm(A - QH*RH, "fro") / norm_A;
    Y_3(m-99) = norm(A - U*D*V, "fro") / norm_A;
end

% LS fitting

lX_1 = log(X_1);
lX_2 = log(X_2);
lX_3 = log(X_3);

lY_1 = log(Y_1);
lY_2 = log(Y_2);
lY_3 = log(Y_3);

p_1 = polyfit(lX_1,lY_1,1);
p_2 = polyfit(lX_2,lY_2,1);
p_3 = polyfit(lX_3,lY_3,1);

fprintf("Least sq soln for MGS: q = %f, c = %f\n", p_1(1), p_1(2));
fprintf("Least sq soln for Householder: q = %f, c = %f\n", p_2(1), p_2(2));
fprintf("Least sq soln for SVD: q = %f, c = %f\n", p_3(1), p_3(2));

y1 = polyval(p_1,lX_1);
figure
plot(lX_1,lY_1,'.');
hold on
plot(lX_1,y1); title('MGS');
hold off

y2 = polyval(p_2,lX_2);
figure
plot(lX_2,lY_2,'.');
hold on
plot(lX_2,y2); title('Householder');
hold off

y3 = polyval(p_3,lX_3);
figure
plot(lX_3,lY_3,'.');
hold on
plot(lX_3,y3); title('SVD');
hold off


% Distribution of Y (as histogram)
Y_1 = zeros(1, 1000); Y_2 = zeros(1, 1000); Y_3 = zeros(1, 1000);
for i=1:3
    m = 100*i;
    for j=1:1000
        A = randn(m,m);
        norm_A = norm(A, "fro");

        [QM,RM]=mgs(A);
        [QH,RH]=qr(A);
        [U, D, V] = svd(A);
        Y_1(j) = norm(A - QM*RM, "fro") / norm_A;
        Y_2(j) = norm(A - QH*RH, "fro") / norm_A;
        Y_3(j) = norm(A - U*D*V, "fro") / norm_A;
    end
    edges = [-10 -0.2:0.05:0.1 10];
    figure; histogram(Y_1, edges); title(['MGS, m=', m, ""]);
    figure; histogram(Y_2, edges); title(['Householder, m=', m, ""]);
    figure; histogram(Y_3); title(['SVD, m=', m, ""]);
end

function [q, r] = mgs(a)
    [m,n] = size(a);
    v = zeros(m,n);
    q = zeros(m,n);
    r = zeros(n,n);
    for i = 1:n
        v(:, i) = a(:, i);
    end
    for i = 1:n
        r(i,i) = norm(v(:, i));
        if r(i,i) ~= 0
            q(:, i) = v(:, i) / r(i,i);
        end
        for j = i+1:n
           r(i,j) = q(:,i)' * v(:,j);
           v(:,j) = v(:,j) - r(i,j)*q(:,i);
        end
    end
end