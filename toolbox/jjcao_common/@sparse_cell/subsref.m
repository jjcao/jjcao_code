% Version: Oct. 2004
% Author : Chen Yanover
% Overloading subref function for the sparse_cell class
% 
% Revised: March 2009
function [varargout] = subsref(a,s)

s1 = s(1);

switch s1.type
 case '()'		
	im = subsref(a.indMat, s1);
	nonEmpty = find(im);
	b.nelem = length(nonEmpty);
	c = a.cells(im(nonEmpty));
	im(nonEmpty) = 1:length(nonEmpty);
	b.indMat = im;
	b.cells = c;
	b = class(b, 'sparse_cell');
	
	varargout = {b};
	
 case '{}'
	s1_ = s1; s1_.type = '()';
	im = subsref(a.indMat, s1_);    
	num_ent = length(im(:));
	
	for i=1:num_ent,
		if ~im(i), varargout{i} = [];
		else,      varargout{i} = a.cells{im(i)};
		end
	end
	
 otherwise
	error;
end

s = s(2:end);
if ~isempty(s),
	b = subsref(varargout{:},s); 
	varargout = {b};
end