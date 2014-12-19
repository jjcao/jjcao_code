function retval = texture_coords_exist(IFSObj,varargin)

retval = ~isempty(IFSObj.UV);
