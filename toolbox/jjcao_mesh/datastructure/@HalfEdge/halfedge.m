function HEObj = halfedge(varargin)

% Copy Constructor:
if length(varargin) == 1 && isa(varargin{1},'HalfEdge')
    HEObj = varargin{1};
    return;
end;

% IndexedFaceSet Constructor:
if length(varargin) == 1 && isa(varargin{1},'IndexedFaceSet')
    IFSObj = varargin{1};
else
    IFSObj = IndexedFaceSet(varargin{:});
end;

% HalfEdge fields:
HEObj.ne  = uint32(0);   % Number of edges
HEObj.EO  = uint32([]);  % Edge orig vertex
HEObj.ED  = uint32([]);  % Edge destination vertext
HEObj.ET  = uint32([]);  % Edge twin
HEObj.EF  = uint32([]);  % Edge face
HEObj.EN  = uint32([]);  % Next edge in the face 
HEObj.EP  = uint32([]);  % Previous edge in the face 

HEObj.VE  = uint32([]);  % One of the vertex orig edges
HEObj.FE  = uint32([]);  % One of the face edges

HEObj.nb  = uint32(0);   % Number of boundaries
HEObj.BFS = uint32([]);  % Boundary faces size (number of vertices in each boundary face)
HEObj.IBF = uint8([]);   % Boundary faces boolean flags

if length(varargin) == 0
    IFSObj = IndexedFaceSet;    
    HEObj = class(HEObj,'HalfEdge',IFSObj);
    return
end

try
    HEObj = p_ifs2he(IFSObj,HEObj);
catch
    error(sprintf('Error constructing HalfEdge object.\n%s', lasterr));
end;
    
HEObj = class(HEObj,'HalfEdge',IFSObj);
