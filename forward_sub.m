function x= forward_sub(L,b)
% Use forward substitution to solve the lower triangular system  L*x=b 

  [m,n]=size(L); 
  x=zeros(size(b));
  for k=1:m
    x(k)=(b(k)-L(k,1:k-1)*x(1:k-1))/L(k,k);
  end
