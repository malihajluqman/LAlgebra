function H=househess(A)

[m, n]=size(A);


for k=1:(m-2)
    x=A(k+1:m,k);
    vk=x;
    vk(1)=sign(x(1))*norm(x)+x(1);
    vk=vk/norm(vk);
    A(k+1:m,k:n)=A(k+1:m, k:n)-2*vk*(vk'*A(k+1:m,k:n));
    A(1:m,k+1:n)=A(1:m,k+1:n)-2*(A(1:m,k+1:n)*vk)*vk';
end
H=A;