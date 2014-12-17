function retval = vertex_dest_edges(HEObj,varargin)

retval = vertex_orig_edges(HEObj,varargin{:});
if ~iscell(retval)
    retval = edge_twin(HEObj,retval);
else
	for i=1:length(retval),
        retval{i} = edge_twin(HEObj,retval{i});
	end;
end;

