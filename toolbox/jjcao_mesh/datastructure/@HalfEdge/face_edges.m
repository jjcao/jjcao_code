function retval = face_edges(HEObj,varargin)

if isempty(varargin) || isempty(varargin{1})
    f = 1:nof(HEObj);
else
    f = varargin{1};
end;

retval = {};

for i=1:length(f),
    cf = f(i);
    oe = [0 HEObj.FE(cf)];
    
    while (oe(end) ~= oe(1))
        oe(1) = oe(2);
        oe = [oe edge_next(HEObj,oe(end))];
    end;
    
    retval{i} = oe(2:end-1);
end;
    
if length(retval) == 1
    retval = retval{1};
end;
