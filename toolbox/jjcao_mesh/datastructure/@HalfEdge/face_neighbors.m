function retval = face_neighbors(HEObj,varargin)

retval = face_edges(HEObj,varargin{:});
if ~iscell(retval)
    retval = edge_face(HEObj,edge_twin(HEObj,retval));
else
	for i=1:length(retval),
        retval{i} = edge_face(HEObj,edge_twin(HEObj,retval{i}));
	end;
end;
