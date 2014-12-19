function retval = face_colors(IFSObj,varargin)

if ~face_colors_exist(IFSObj)
    error('Face colors do not exist.');
end;

if isempty(varargin) || isempty(varargin{1})
    f = 1:nof(IFSObj);
else
    f = varargin{1};
end;

retval = double(IFSObj.FC(f,:))/255;
