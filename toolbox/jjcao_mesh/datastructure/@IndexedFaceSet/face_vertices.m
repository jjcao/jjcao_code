function retval = face_vertices(IFSObj,varargin)

if isempty(varargin) || isempty(varargin{1})
    f = 1:nof(IFSObj);
else
    f = varargin{1};
end;

retval = {IFSObj.F{f}};
for i=1:length(retval),
    retval{i} = double(retval{i});
end;

if length(retval) == 1
    retval = retval{1};
end;


