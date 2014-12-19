function retval = face_normals(IFSObj,varargin)

if ~face_normals_exist(IFSObj)
    error('Face normals do not exist.');
end;

if isempty(varargin) || isempty(varargin{1})
    f = 1:nof(IFSObj);
else
    f = varargin{1};
end;

retval = double(IFSObj.FN(f,:));
