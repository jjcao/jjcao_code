function retval = boundary_edges(HEObj,varargin)

if isempty(varargin) || isempty(varargin{1})
    b = 1:nob(HEObj);
else
    b = varargin{1};
end;

retval = {};
for i = 1:length(b),
    cbe = double(HEObj.FE(nof(HEObj.indexedfaceset)+b(i)));
    nbe = edge_next(HEObj,cbe);
    while cbe ~= nbe
        cbe = [cbe nbe];
        nbe = edge_next(HEObj,cbe(end));
    end;
    retval{i} = cbe;
end;

if length(retval) == 1
    retval = retval{1};
end;


