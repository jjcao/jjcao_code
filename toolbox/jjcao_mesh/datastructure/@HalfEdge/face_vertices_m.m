function retval = face_vertices_m(HEObj,varargin)

if isempty(varargin) || isempty(varargin{1})
    ifci = 1:nof(HEObj.indexedfaceset);
    bfci = 1:nob(HEObj);
else
    f = varargin{1};
    ifci = f(find(f<=nof(HEObj.indexedfaceset)));
    bfci = setdiff(f,ifci);  bfci = bfci - nof(HEObj.indexedfaceset);
end;

ifc = face_vertices_m(HEObj.indexedfaceset,ifci);
bfc = boundary_vertices_m(HEObj,bfci);

if length(ifci) == 0
    retval = bfc;
    return;
end;

if length(bfci) == 0
    retval = ifc;
    return;
end;

[in,im] = size(ifc);
[bn,bm] = size(bfc);

retval = zeros(in+bn,max(im,bm));

if im == bm
    retval(1:in,:) = ifc;
    retval(in+1:end,:) = bfc;
elseif im < bm
    retval(1:in,1:im) = ifc;
    retval(1:in,im+1:end) = ifc(:,1)*ones(1,bm-im);
    retval(in+1:end,:) = bfc;
else
    retval(1:in,:) = ifc;
    retval(in+1:end,1:bm) = bfc;
    retval(in+1:end,bm+1:end) = bfc(:,1)*ones(1,im-bm);
end;
