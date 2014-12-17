function retval = is_boundary_face(HEObj,varargin)

if isempty(varargin) || isempty(varargin{1})
    f = 1:nof(HEObj);
else
    f = varargin{1};
end;

retval = double(HEObj.IBF(f)) == 1;
