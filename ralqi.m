%Computes one eigenvalue  and eigenvector
%of A using the Rayleigh quotient
%iteration, close to initial eigenvector estimate

function [s,v, k]=ralqi(A,v0, tol)
format long

v=v0/norm(v0);
s=v'*A*v;
[m,n]=size(A);
sk=s+1;
k=1;

while abs(sk-s)>tol
    [P, L, U]=gausspivot(A-s*eye(m));
    w=back_sub(U,forward_sub(L,P*v));
    v=w/norm(w);
    sk=s;
    s=v'*A*v;
    k=k+1;
end