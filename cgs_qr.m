function [Q,R] = cgs_qr(A)
% Computation of the "skinny" QR decomposition of an m-by-n matrix
% (m>=n) using the the Classical Gram-Schmidt algorithm.
% See Trefethen and Bau, Algorithm 8.1
[m,n] = size(A);
if m<n, fprintf('ERROR: A should be an m-by-n matrix with m >= n.\n'); end
Q = zeros(m,n);
R = zeros(n,n);
for k=1:n
Q(:,k) = A(:,k);
for j=1:k-1
R(j,k) = Q(:,j)'*A(:,k);
Q(:,k) = Q(:,k) - R(j,k)*Q(:,j);
end
R(k,k) = norm(Q(:,k));
Q(:,k) = Q(:,k)/R(k,k);
end
