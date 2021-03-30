## The function fast_shor_algorithm_procedure implements the quantum part 
## of the algorithm described in 3-7 of the Home Assignment 4. It takes an 
## integer smaller than 256 as input, and outputs a guess of s.


function s = fast_shor_algorithm_procedure(m);
%% initialize 21-qubit state q:
  % initialize the 13-qubit state in superposition, superpos_work_qubits.
  work_qubits = ket('0000000000000');
  superpos_work_qubits = had(work_qubits);
  
  % initialize 8-qubit, zero state, helper_qubits  
  helper_qubits = ket('00000000');
  
  % Build the matrix representation of state q, called Q
  tranpose_helper_qubits = helper_qubits.'; 
  Q = superpos_work_qubits*tranpose_helper_qubits;
  
%% Power function applied to Q: 
  % Choose remainder r in multiplicative group Z_m  
  r = multiplicative_remainder(m);
  
  % Apply power function with base r
  d = 2^13;
  r_a = 1;
  for a = 0:d-1;
    r_a = mod(r_a*r, m);
    Q(a+1, r_a+1) = Q(a+1, 1);
    Q(a+1, 1) = 0;
  endfor

%% 8-bit measurement u: 
  % sum each column of Q and hold results in 2^13-dimensional vector, v
  dimv = 2^8;
  v = zeros(1,dimv);
  for j =1:dimv;
    v(j) = sum(Q(1:2^13, j));
  endfor

  
  % Normalize vector of sums, and produce vector of measurement probabilities v
  mag_v = sqrt(sumsq(v));
  v = v/mag_v;
  v = v.*v;
  
  % Monte Carlo method using cummulative sum of measurement probabilities in v
  t8 = rand;
  vsum = cumsum(v);
  u = 0;
  i = 1;
  while i < dimv+1;
    if t8 < vsum(i);
      u = i; 
      i = dimv+1;
    endif
    i = i +1;
  endwhile
## It is worth noting that although we set u = i in the above, the actual 
## measurement outcome is i-1. Octave is 1-indexed so some adjustments are 
## made for the final output, however the program works with u = i.
  
%% Fourier transform the column of Q corresponding to the measurement u:
  qfft = fft(Q(1:2^13, u));

%% 13-bit measurement c (similar 8-bit measurement):
  % Normalize qfft, and produce a vector of measurement probabilities qfft_prob 
  qfft_conj = conj(qfft);
  qfft_norm = qfft.*qfft_conj;
  mag_qfft = sqrt(sum(qfft_norm));  
  qfft_norm = sqrt(qfft_norm)/mag_qfft;
  qfft_prob = qfft_norm.*qfft_norm; 
   
  % Monte Carlo method using cummulutive sum of qfft_prob 
  qsum = cumsum(qfft_prob);
  dimq = length(qsum);
  t13 = rand;
  c = 0;
  l = 1;
  while l < dimq+1;
    if t13 < qsum(l);
      c = l;
      l = dimq+1;
    endif
    l = l+1;
  endwhile
  
%% Compute distance from zc/2^13 to nearest int; guess s:
  z_min = ceil(m-3*sqrt(m));
  z_max = floor(m+1-2*sqrt(m));
  dist = 1;
  s = 0;
  c = c-1; % adjusting for 1-indexing of octave 
  for z = z_min:z_max;
    dist_temp = abs(z*c/(2^13) - round(z*c/(2^13)));
    if dist_temp < dist;
      dist = dist_temp;
      s = z;
    endif
  endfor
endfunction