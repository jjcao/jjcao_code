function retval = is_boundary_edge(HEObj,varargin)

retval = is_boundary_face(HEObj,edge_face(HEObj,varargin{:}));
