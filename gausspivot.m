function [P, L, U]=gausspivot(A)

[m, n]=size(A);
if m~=n disp('The matrix must be an mxm matrix')
end


U=A;
L=eye(m);
P=eye(m);

for k=1:(m-1)
    M=0;
    %Finding the maximum value in the column k
    for i=k:m
        if abs(U(i,k))>M
            M=abs(U(i,k));
        end
    end
    %Switching rows 
    if M > U(k,k)
        tempu=U(k,k:m);
        U(k,k:m)=U(i,k:m);
        U(i,k:m)=tempu;
        templ=L(k,1:(k-1));
        L(k,1:(k-1))=L(i,1:(k-1));
        L(i,1:(k-1))=templ;
        tempp=P(k,:);
        P(k,:)=P(i,:);
        P(i,:)=tempp;    
    end
    %Performing Gaussian elimination
    for j=(k+1):m
        L(j,k)=U(j,k)/U(k,k);
        U(j,k:m)=U(j,k:m)-L(j,k)*U(k,k:m);
    end
end
    
    
        