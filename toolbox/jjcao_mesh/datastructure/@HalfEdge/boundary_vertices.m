function retval = boundary_vertices(HEObj,varargin)

retval = boundary_edges(HEObj,varargin{:});
if ~iscell(retval)
    retval = edge_orig(HEObj,retval);
else
	for i=1:length(retval),
        retval{i} = edge_orig(HEObj,retval{i});
	end;
end;


