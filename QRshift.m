function [V, L, k]=QRshift(A, tol)

format long
H=househess(A);
A=H;
[m,n]=size(A);
mu=A(m,m);
[Q, R]=mgs_qr(A-mu*eye(m));
A=R*Q+mu*eye(m);
V=Q;
k=1;

while norm(A-diag(diag(A)))>tol
    k=k+1;
    mu=A(m,m);
    [Q, R]=mgs_qr(A-mu*eye(m));
    A=R*Q+mu*eye(m);
    V=V*Q;
end
L=diag(diag(A));

    