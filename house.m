function [W,R] = house(A)
% Compute the QR factorization of real A using Householder reflectors
% See Trefethen and Bau, Numerical Linear Algebra, Algorithm 10.1 (page 73)
% The object V is a "cell array": vk (a vector of length m-k+1) is stored in V{k}.
[m,n] = size(A);
sz=size(A);
W=zeros(sz);
for k=1:min(m-1,n)
ak = A(k:m,k); % vector to be zeroed out
vk = ak; vk(1) = vk(1) + sign(ak(1))*norm(ak); % 1. construct vk that defines the reflector
vk = vk/norm(vk); % 2. normalize vk
A(k:m,k:n) = A(k:m,k:n) - 2*vk*(vk'*A(k:m,k:n)); % 3. update A
W(k:m,k)=vk;
end

R = A;