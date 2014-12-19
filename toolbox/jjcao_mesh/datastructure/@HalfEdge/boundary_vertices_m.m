function retval = boundary_vertices_m(HEObj,varargin)

if isempty(varargin) || isempty(varargin{1})
    b = 1:nob(HEObj);
else
    b = varargin{1};
end;

bv = boundary_vertices(HEObj,b);
if ~iscell(bv)
    bv = {bv};
end;
maxfs = max(double(HEObj.BFS(b)));
retval = zeros(length(b),maxfs);

for i=1:length(b)
    retval(i,:) = [bv{i} bv{i}(1)*ones(1,maxfs-double(HEObj.BFS(b(i))))];
end;

