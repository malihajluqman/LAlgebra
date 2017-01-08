function [Q,R] = mgs_qr(A)
% Computation of the "skinny" QR decomposition of an m-by-n matrix
% (m>=n) using the the Modified Gram-Schmidt algorithm.
% See Trefethen and Bau, Algorithm 8.1
[m,n] = size(A);
if m<n, fprintf('ERROR: A should be an m-by-n matrix with m >= n.\n'); end
Q = zeros(m,n);
R = zeros(n,n);
Q = A;
for j=1:n
R(j,j) = norm(Q(:,j));
Q(:,j) = Q(:,j)/R(j,j);
for k=j+1:n
R(j,k) = Q(:,j)'*Q(:,k);
Q(:,k) = Q(:,k) - R(j,k)*Q(:,j);
end
end
