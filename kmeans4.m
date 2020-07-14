function [C] = kmeans4(I)
k=4;
[m, n] = size(I);
X = reshape(double(I), m*n, 1);
A = X;
A(find(A==0))=[];
[a, b] = size(A);
rng('default');
C = A(randperm(a*b, k), :);
J_prev = inf; iter = 0; J = []; tol = 1e-2;
while true,
    iter = iter + 1;
    dist = sum(X.^2, 2)*ones(1, k) + (sum(C.^2, 2)*ones(1, m*n))' - 2*X*C';
    [~, label] = min(dist, [], 2) ;
    for i = 1:k,
       C(i, :) = mean(X(label == i , :));
    end
    J_cur = sum(sum((X - C(label, :)).^2, 2));
    J = [J, J_cur];
    sprintf('#iteration: %03d, objective fcn: %f', iter, J_cur);
    if norm(J_cur-J_prev, 'fro') < tol,
        break;
    end
    J_prev = J_cur;
end
C=sort(C);

