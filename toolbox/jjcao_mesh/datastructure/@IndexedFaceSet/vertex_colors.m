function retval = vertex_colors(IFSObj,varargin)

if ~vertex_colors_exist(IFSObj)
    error('Vertex colors do not exist.');
end;

if isempty(varargin) || isempty(varargin{1})
    v = 1:nov(IFSObj);
else
    v = varargin{1};
end;

retval = double(IFSObj.VC(v,:))/255;
