function R=cholesky(A)
R=A;
[m, n]=size(A);
if m~=n disp('Matrix must be a square matrix')
end

for k=1:m
    for j=(k+1):m
        R(j,j:m)=R(j,j:m)-R(k,j:m)*conj(R(k,j))/R(k,k);
    end
    R(k,k:m)=R(k,k:m)/sqrt(R(k,k));
    for i=1:(k-1)
        R(k,i)=0;
    end
end