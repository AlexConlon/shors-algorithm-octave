function p = ket(s);
  l = length(s); 
  p = 1;
  for i = 1:l;
    c = comp_vect(bin2dec(s(i)));
    p = kron(p, c); 
    i = i +1;
  endfor
endfunction
