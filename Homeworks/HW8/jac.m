
function [Df x] = jac (x, f) 

  n = numel(x);
	 
  h=1e-6;
  
  fixed = false(n,1);

  x = repmat (x(:), 1, n) + h * 1i * eye (n);

  idx = find (! fixed);

  ## after first evaluation, dimensionness of 'f' is known
  t_Df = imag (f (x(:, idx(1)))(:));
  dim = numel (t_Df);

  Df = zeros (dim, n);

  Df(:, idx(1)) = t_Df;

  for count = idx(2:end)'
    Df(:, count) = imag (f (x(:, count))(:));
  endfor

  Df /=  h;

endfunction
