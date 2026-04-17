function dX = deriv1_cols(X, dt)
[n,K] = size(X);
dX = zeros(n,K);

if K < 2
    return;
end

dX(:,1) = (X(:,2) - X(:,1)) / dt;
dX(:,K) = (X(:,K) - X(:,K-1)) / dt;

for k = 2:K-1
    dX(:,k) = (X(:,k+1) - X(:,k-1)) / (2*dt);
end
end
