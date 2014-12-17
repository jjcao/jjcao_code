function retval = face_size(IFSObj,varargin)

if isempty(varargin) || isempty(varargin{1})
    f = 1:nof(IFSObj);
else
    f = varargin{1};
end;

retval = double(IFSObj.FS(f));
