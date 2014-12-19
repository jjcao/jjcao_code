% Version: Oct. 2004
% Author : Chen Yanover
% Overloading length function for the sparse_cell class
% 
% Revised: March 2009
function k = length(sc)

k = length(sc.indMat);