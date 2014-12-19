% Version: March 2009
% Author : Chen Yanover
% Overloading numel function for the sparse_cell class
function n = numel(a,varargin)

n = numel(a.indMat, varargin{:});