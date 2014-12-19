function retval = vertex_orig_edges(HEObj,varargin)

if isempty(varargin) || isempty(varargin{1})
    v = 1:nov(HEObj);
else
    v = varargin{1};
end;

retval = {};

for i=1:length(v),
    cv = v(i);
    oe = [0 HEObj.VE(cv)];
    
    while (oe(end) ~= oe(1))
        oe(1) = oe(2);
        oe = [oe edge_next(HEObj,edge_twin(HEObj,oe(end)))];
    end;
    
    retval{i} = double(oe(end:-1:3));
end;
    
if length(retval) == 1
    retval = retval{1};
end;
