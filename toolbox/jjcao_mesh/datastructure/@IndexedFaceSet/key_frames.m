function retval = key_frames(IFSObj,varargin)

if ~key_frames_exist(IFSObj)
    error('Key frames do not exist.');
end;

if isempty(varargin) || isempty(varargin{1})
    a = 1:length(IFSObj.A);
else
    a = varargin{1};
end;

retval = {IFSObj.A{a}};
for i=1:length(retval),
    retval{i} = double(retval{i});
end;

if length(retval) == 1
    retval = retval{1};
end;


