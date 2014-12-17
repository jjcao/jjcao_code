function retval = face_vertices_m(IFSObj,varargin)

if isempty(varargin) || isempty(varargin{1})
    f = 1:nof(IFSObj);
else
    f = varargin{1};
end;

maxvf = max(face_size(IFSObj,f));
retval = zeros(length(f),maxvf);

for i=1:length(f)
    retval(i,:) = [double(IFSObj.F{f(i)}) double(IFSObj.F{f(i)}(1))*ones(1,maxvf-face_size(IFSObj,f(i)))];
end;

