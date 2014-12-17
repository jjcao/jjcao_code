% Version: Nov. 2004
% Author : Chen Yanover
% Overloading subassign function for the sparse_cell class
% 
% Revised: March 2009
function a = subsasgn(a,s,b)

s1 = s(1);
s  = s(2:end);

[m,n] = size(a);

n_subs = length(s1.subs);
if n_subs>2, error('ND sparce cell is not supported'); end

switch n_subs, 
 case 1, 
	num_cells = prod(size(a));
	if strcmp(s1.subs{1},':'), 
		s1.subs{1} = 1:num_cells;
	elseif max(s1.subs{1})>num_cells,
		a.indMat = a.indMat(:)';
		a.indMat(max(s1.subs{1})) = 0;
	end
	subs = s1.subs{1};
	
 case 2, 
	for i=1:n_subs, 	 
		if strcmp(s1.subs{i},':'), 
			s1.subs{i} = 1:size(a.indMat,i);
		end
		max_ind(i) = max([s1.subs{i}, size(a.indMat,i)]);
	end

	if max_ind(1)>m | max_ind(2)>n,
		a.indMat(max_ind(1), max_ind(2)) = 0;
		[m,n] = size(a.indMat);
	end
	
	% trasnform subrefs to inds
	subs1 = s1.subs{1}(:)*ones(1,length(s1.subs{2}));
	subs2 = (s1.subs{2}(:)*ones(1,length(s1.subs{1})))';
	subs  = sub2ind([m, n], subs1(:), subs2(:)); 
end

switch s1.type, 
 case '{}', 	

	if length(subs)>1, 
		error('Cell assignment to multiple entries is not supported');
	end
		
	if isempty(s),
		a = update_sinle_cell(a, subs, b);		
	else,
		ind = a.indMat(subs);
		a.cells{ind} = subsasgn(a.cells{ind}, s, b);
	end		
	
 case '()',	
	if length(b)~=1 & length(b)~=length(subs), 
		error('Subscripted assignment dimension mismatch');
	end
	
	for i=1:length(subs), % now 1d ...
		
		b_ind = i;
		if length(b)==1, b_ind=1; end
		
		if isa(b, 'sparse_cell'),

			b_i = subsref(b, substruct('()', {b_ind}));
			bs  = substruct('{}', {':'});
			a   = update_sinle_cell(a, subs(i), subsref(b_i,bs));			
		else,
			b_i = b(b_ind);
			a = update_sinle_cell(a, subs(i), b_i{:});		
		end
	end

 otherwise,
	error('sparse_cell, subasgn: wrong usage');
end


function [a,ind] = update_sinle_cell(a, i, b)

ind = a.indMat(i);
if ~ind & isempty(b),       % nothing to do ... 
	;
	
elseif ind  &  isempty(b),  % remove the existing element 
	a.indMat(i)          = 0;
	a.cells(ind)         = [];
	a.nelem              = a.nelem-1;
	l_inds               = find(a.indMat>ind);
	a.indMat(l_inds) = a.indMat(l_inds) - 1; 
	
elseif ~ind & ~isempty(b),  % add a cell
	a.indMat(i)      = a.nelem+1;			
	a.cells{end+1}   = b;
	a.nelem          = a.nelem + 1;
	
elseif ind  & ~isempty(b),  % just update the cell 
	a.cells{ind}   = b;
end
