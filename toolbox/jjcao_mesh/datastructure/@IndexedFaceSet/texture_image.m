function retval = texture_image(IFSObj,varargin)

if ~texture_image_exist(IFSObj)
    error('Texture image do not exist.');
end;
    
retval = IFSObj.TImg;
