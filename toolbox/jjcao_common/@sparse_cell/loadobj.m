function a = loadobj(a)

if isstruct(a), 
	a = rmfield(a, {'n','m'});
	a = class(a, 'sparse_cell');
	
	a.nelem = length(a.cells);
	if sum(a.indMat(:)~=0)~=a.nelem,
		error('load failed');
	end
end