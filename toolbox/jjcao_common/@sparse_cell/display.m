% Version: Oct. 2004
% Author : Chen Yanover
% Overloading display function for the sparse_cell class
% 
% Revised: March, 2009
function display(sc)

disp(' ');
disp([inputname(1) ' = ']);
disp(' ');

if ~sc.nelem, 
	fprintf('   All empty sparse_cell: %d-by-%d\n\n', size(sc));
else
	[mi,ni] = find(sc.indMat');

	for i=1:length(ni),
		ind_str = sprintf('(%d,%d) ', ni(i), mi(i));
		fprintf('%11s  ', ind_str);
		
		cell_i  = sc.cells{sc.indMat(ni(i),mi(i))};				
		
		str_i = [];
		if length(cell_i)>10 || all(size(cell_i)>1),
			size_i = size(cell_i);						
			
			for d=1:ndims(cell_i)-1,
				str_i = [str_i, sprintf('%dx', size_i(d))];
			end
			str_i = [str_i, sprintf('%d %5s', size_i(end), class(cell_i))]; 			
			fprintf('[%12s]\n', str_i);
		else
			switch class(cell_i),
			 case 'double',
				str_i = [str_i, sprintf(' %.3f', cell_i)];
				fprintf('[%12s ]\n', str_i);
			 case {'uint8', 'uint16', 'uint32', 'uint64', 'int8', 'int16', 'int32', 'int64'}
				str_i = [str_i, sprintf(' %3d', cell_i)];
				fprintf('[%12s ]\n', str_i);
			 case 'char', 
				str_i = [str_i, sprintf('%s', cell_i)];
				fprintf('[%12s]\n', str_i);
			 otherwise, 
				disp(cell_i);
			end
		end		

	end
end

