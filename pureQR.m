function [V, L, k]=pureQR(A, tol)
format long
[Q, R]=mgs_qr(A);
A=R*Q;
V=Q;
k=1;

while norm(A-diag(diag(A)))>tol
    k=k+1;
    [Q, R]=mgs_qr(A);
    A=R*Q;
    V=V*Q;
end
L=diag(diag(A));