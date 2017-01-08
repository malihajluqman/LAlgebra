function x=back_sub(U,b)
% Use backward substitution to solve the upper triangular system  U*x=b 

  [m,n]=size(U); 
  x=zeros(size(b));
  for k=m:-1:1
    x(k)=(b(k)-U(k,k+1:m)*x(k+1:n))/U(k,k);
  end
