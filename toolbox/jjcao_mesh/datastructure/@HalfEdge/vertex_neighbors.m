function retval = vertex_neighbors(HEObj,varargin)

retval = vertex_orig_edges(HEObj,varargin{:});
if ~iscell(retval)
    retval = edge_dest(HEObj,retval);
else
	for i=1:length(retval),
        if ~isempty(retval{i})
            retval{i} = edge_dest(HEObj,retval{i});
        end
	end
end
