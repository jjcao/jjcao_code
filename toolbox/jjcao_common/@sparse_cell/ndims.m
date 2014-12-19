% Version: March 2009
% Author : Chen Yanover
% Overloading ndims function for the sparse_cell class
function n = ndims(a)

n = ndims(a.indMat);