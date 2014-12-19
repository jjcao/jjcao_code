function retval = face_size(HEObj,varargin)

if isempty(varargin) || isempty(varargin{1})
    ifci = 1:nof(HEObj.indexedfaceset);
    bfci = 1:nob(HEObj);
else
    f = varargin{1};
    ifci = f(find(f<=nof(HEObj.indexedfaceset)));
    bfci = setdiff(f,ifci);  bfci = bfci - nof(HEObj.indexedfaceset);
end;

retval = [];
if ~isempty(ifci)
    retval = [retval face_size(HEObj.indexedfaceset,ifci)];
end;

retval = [retval double(HEObj.BFS(bfci));];

