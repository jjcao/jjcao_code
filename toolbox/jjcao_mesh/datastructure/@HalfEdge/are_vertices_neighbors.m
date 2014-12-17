function retval = are_vertices_neighbors(HEObj,VL,varargin)

if ~isnumeric(VL) || length(VL) < 2
    error('Vertex-List is expected to be a vector with at least two vertex indices.');
end;

retval = 0;
for i=1:length(VL)-1,
    vn = vertex_neighbors(HEObj,VL(i));
    if length(intersect(vn,VL(i+1:end))) ~= length(VL)-i
        return;
    end;
end;
retval = 1;