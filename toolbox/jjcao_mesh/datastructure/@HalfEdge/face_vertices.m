function retval = face_vertices(HEObj,varargin)

if isempty(varargin) || isempty(varargin{1})
    ifci = 1:nof(HEObj.indexedfaceset);
    bfci = 1:nob(HEObj);
else
    f = varargin{1};
    ifci = f(find(f<=nof(HEObj.indexedfaceset)));
    bfci = setdiff(f,ifci);  bfci = bfci - nof(HEObj.indexedfaceset);
end;

if ~isempty(ifci)
    ifc = face_vertices(HEObj.indexedfaceset,ifci);
	if ~iscell(ifc)
        ifc = {ifc};
	end;
else
    ifc = {};
end;

if ~isempty(bfci)
    bfc = boundary_vertices(HEObj,bfci);
	if ~iscell(bfc)
        bfc = {bfc};
	end;
else
    bfc = {};
end;

retval = {ifc{:} bfc{:}};

if length(retval) == 1
    retval = retval{1};
end;
