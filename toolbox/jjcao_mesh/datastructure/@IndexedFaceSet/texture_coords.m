function retval = texture_coords(IFSObj,varargin)

if ~texture_coords_exist(IFSObj)
    error('Texture coordinates do not exist.');
end;
    
if isempty(varargin) || isempty(varargin{1})
    v = 1:nov(IFSObj);
else
    v = varargin{1};
end;

retval = double(IFSObj.UV(v,:));
