## As inputs, the function Shors_algorithm takes an integer m (with only two 
## prime factors) to factor, and the number of runs, n, of 
## fast_shor_algorithm_procedure. For n = 100 it takes about 100 seconds to run 
## this program it beings by initialize an n-dimensional vector with all zero 
## entries called s_array. This vector will hold the results of each run.

function factors = Shors_algorithm(m, n);   
  s_array = zeros(1,n);
  
%% for-loop used to populate s_array with a guess of s for each run 
  for i = 1:n;
    s_array(i) = fast_shor_algorithm_procedure(m);
  endfor

## mode() returns the element of s_array with the most occurrences
  s = mode(s_array);
  
## finding roots of x^2-(m+1-s)+m
  p_1 = ((m+1-s)+sqrt((m+1-s)^2-4*m))/2;
  p_2 = ((m+1-s)-sqrt((m+1-s)^2-4*m))/2;
  
## returns (in theory) the prime factors of m
  factors = [p_1, p_2];
endfunction
