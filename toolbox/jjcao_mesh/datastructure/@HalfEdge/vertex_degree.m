function retval = vertex_degree(HEObj,varargin)

vn = vertex_neighbors(HEObj,varargin{:});
if ~iscell(vn)
    retval = length(vn);
else
    retval = zeros(1,length(vn));
    for i=1:length(vn)
        retval(i) = length(vn{i});
    end;
end;

