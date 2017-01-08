%Computes one eigenvalue of A using power iteration
%Outputs: s, v are the computed eigenvalue and eigenvector
%         that is close to the initial eigenvector estimate
%         k is the number of iterations 
%Inputs: v0 is the initial estimate of the eigenvector
%        tol is the user specified tolerance

function [s, v, k]=poweriter(A,v0, tol)
format long

k=1;
w=A*v0;
v=w/norm(w);
s=v'*A*v;
sk=s+1;

while abs(sk-s)>tol
    k=k+1;
    w=A*v;
    v=w/norm(w);
    sk=s;
    s=v'*A*v;
end