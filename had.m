function qx = had(q);
  # I had to build my own Hadamard function as Octave's Hadamard function didn't
  # support the size I needed.
  l = log2(length(q));
  H0 = (1/sqrt(2))*hadamard(2);
  H = H0;
  
  for i = 1:l-1;
    H = kron(H, H0);
  endfor
  qx = H*q;
endfunction
