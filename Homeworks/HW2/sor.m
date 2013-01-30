function x = sor(a, w, b, x)

for i=1:10
  n = length(a);
  k = size(x,2)+1;
  x(:,k) = zeros(n,1);

#  for i = 1:n
#      
#    x(i,k) = (1-w)*x(i,k-1) + (w/a(i,i))*(b(i) -
#					  sum(a(i,1:i).*x(:,k)) -
#					  sum(a(i,i+1:n).*x(:,k-1)));
#  endfor

  D = eye(length(a)).*a;
  L = tril(a, -1);
  U = triu(a,  1);

  T = inverse(D-w*L)*((1-w)*D + w*U);
  c = w*inverse(D-w*L)*b;

#  x(:, k) = inverse(D+w*L)*(w*b - (w*U + (w-1)*D)*x(:,k-1))

  x(:,k) = T*x(:,k-1) + c;

endfor
endfunction

