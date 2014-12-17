function retval = vertex_normals(IFSObj,varargin)

if ~vertex_normals_exist(IFSObj)
    error('Vertex normals do not exist.');
end;

if isempty(varargin) || isempty(varargin{1})
    v = 1:nov(IFSObj);
else
    v = varargin{1};
end;

retval = double(IFSObj.VN(v,:));
