% Version: Oct. 2004
% Author : Chen Yanover
% Overloading size function for the sparse_cell class 
% 
% Modified by Hoi Wong (Mar 9, 2009)
function varargout = size(sc,varargin)

mat_size = size(sc.indMat, varargin{:});

if nargout<2, 
	varargout = {mat_size};
elseif nargout==2, 
	varargout = mat2cell(mat_size, 1, [1,1]);
end