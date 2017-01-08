function [L, U]=gauss(A)

[m, n]=size(A);
if m~=n disp('The matrix must be an mxm matrix') 
end
U=A;
L=eye(m);
for j=1:m-1
    for i=(j+1):m
        L(i,j)=U(i,j)/U(j,j);
        U(i,j:m)=U(i,j:m)-L(i,j)*U(j,j:m);
    end
end