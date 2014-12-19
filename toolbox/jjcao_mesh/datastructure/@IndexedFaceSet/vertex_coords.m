function retval = vertex_coords(IFSObj,varargin)

if isempty(varargin) || isempty(varargin{1})
    v = 1:nov(IFSObj);
else
    v = varargin{1};
end;

retval = double(IFSObj.V(v,:));
