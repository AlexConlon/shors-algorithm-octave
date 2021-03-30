function t = multiplicative_remainder(m);
% returns an element of the multiplicative group modulo m
  % Chooses a t between 0 and m-1, at random
  t = round(rand*(m-1));
  
  % Verifies if t belongs to the multiplicative group modulo m
  if gcd(t,m) > 1;
    t = multiplicative_remainder(m);
  elseif t == 0;
    t = multiplicative_remainder(m);
  endif
endfunction
