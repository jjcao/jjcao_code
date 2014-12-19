function IFSObj = indexedfaceset(varargin)

% Copy Constructor:
if length(varargin) == 1 && isa(varargin{1},'IndexedFaceSet')
    IFSObj = varargin{1};
    return;
end;

% IndexedFaceSet fields:
IFSObj.nv   = uint32(0);  % Number of vertices
IFSObj.nf   = uint32(0);  % Number of faces
IFSObj.V    = single([]); % Vertices
IFSObj.F    = uint32([]); % Faces
IFSObj.FS   = uint32([]); % Face size (number of vertices in each face)
IFSObj.VN   = single([]); % Vertices normal
IFSObj.FN   = single([]); % Faces normal
IFSObj.VC   = uint8([]);  % Vertices color
IFSObj.FC   = uint8([]);  % Faces color
IFSObj.UV   = double([]); % Texture coordinates
IFSObj.TImg = uint8([]);  % Texture image
IFSObj.A    = cell(0);    % Animation key-frames

if length(varargin) == 0
    IFSObj = class(IFSObj,'IndexedFaceSet');
    return
end

% Parse the constructor command arguments:
if length(varargin) < 2
    error('IndexedFaceSet geometry and connectivity are mandatory.');
end;

% Handling the case of: IndexedFaceSet(V,F)
istart = 1;
if isnumeric(varargin{1}) && (isnumeric(varargin{2}) || iscell(varargin{2}))
    varargin{length(varargin)+1} = 'vertices';
    varargin{length(varargin)+1} = varargin{1};
    varargin{length(varargin)+1} = 'faces';
    varargin{length(varargin)+1} = varargin{2};
    istart = 3;
end;

for i = istart:2:length(varargin)
    if ~ischar(varargin{i})
        error('Parameter name is expected.');
    end;
    
    switch lower(varargin{i})
        
    % Handle the vertices-list/geometry:
    case {'vert','verti','vertic','vertice','vertices','geo','geom','geome','geomet','geometr','geometry'}
        if isempty(varargin{i+1}) || ~isnumeric(varargin{i+1})
            error('Vertices/Geometry value must be supplied as a matrix with vertices position.');
        end;
        IFSObj.V = varargin{i+1};
        
        if size(IFSObj.V,2) > 3
            error('Maximum vertex dimension is 3.');
        end;
        
        if size(IFSObj.V,2) < 3
            IFSObj.V = [IFSObj.V zeros(size(IFSObj.V,1),3-size(IFSObj.V,2))];
        end;
        
        IFSObj.V = single(IFSObj.V);
        IFSObj.nv = size(IFSObj.V,1);
        
    % Handle the faces-list/connectivity
    case {'face','faces','conn','conne','connec','connect','connecti','connectiv','connectivi','connectivit','connectivity'}
        if isempty(varargin{i+1}) 
            error('Faces/Connectivity value must be supplied as a cell array or a matrix.');
        end;
        IFSObj.F = varargin{i+1};
        
        if isa(IFSObj.F,'cell')
            for i=1:length(IFSObj.F)
                IFSObj.FS(i) = uint32(length(IFSObj.F{i}));
                IFSObj.F{i} = uint32(IFSObj.F{i});
            end;
        else                
            f = IFSObj.F;
            IFSObj.F = cell(1,size(f,1));
            for i=1:size(f,1),
                IFSObj.FS(i) = uint32(length(f(i,:)));
                IFSObj.F{i} = uint32(f(i,:));
            end;
            clear f
        end;
        
        IFSObj.nf = length(IFSObj.F);      
        clear i
        
    % Hanlde vertices normals
    case {'vnorm','vnorma','vnormal','vnormals','vertexnormal','vertexnormals'}
        if isempty(varargin{i+1}) || ~isnumeric(varargin{i+1})
            error('Vertex normal value must be supplied as a matrix.');
        end;        
        IFSObj.VN = varargin{i+1};
        
        if size(IFSObj.VN,2) > 3
            error('Maximum normal dimension is 3.');
        end;
        
        if size(IFSObj.VN,2) < 3
            IFSObj.VN = [IFSObj.VN zeros(size(IFSObj.VN,1),3-size(IFSObj.VN,2))];
        end;
        
		nl = sqrt(sum((IFSObj.VN.^2)')'); % Normal vectors length
		nl(nl < 1e-6) = 1; % Not normalize zero length vectors
        IFSObj.VN = IFSObj.VN./(nl*[1 1 1]); % Normalize
        IFSObj.VN = single(IFSObj.VN);
        clear nl
		
    % Hanlde faces normals
    case {'fnorm','fnorma','fnormal','fnormals','facenormal','facenormals'}
        if isempty(varargin{i+1}) || ~isnumeric(varargin{i+1})
            error('Face normal value must be supplied as a matrix.');
        end;        
        IFSObj.FN = varargin{i+1};
        
        if size(IFSObj.FN,2) > 3
            error('Maximum normal dimension is 3.');
        end;
        
        if size(IFSObj.FN,2) < 3
            IFSObj.FN = [IFSObj.FN zeros(size(IFSObj.FN,1),3-size(IFSObj.FN,2))];
        end;
        
		nl = sqrt(sum((IFSObj.FN.^2)')'); % Normal vectors length
		nl(nl < 1e-6) = 1; % Not normalize zero length vectors
        IFSObj.FN = IFSObj.FN./(nl*[1 1 1]); % Normalize
        IFSObj.FN = single(IFSObj.FN);

    % Handle vertices colors:
    case {'vcol','vcolo','vcolor','vcolors','vertexcolor','vertexcolors'}
        if isempty(varargin{i+1}) || ~isnumeric(varargin{i+1})
            error('Vertex color value must be supplied as a matrix.');
        end;        
        IFSObj.VC = varargin{i+1};

        if ~isa(IFSObj.VC,'uint8')
            IFSObj.VC = uint8(255*IFSObj.VC);
        end;

        if size(IFSObj.VC,2) == 3
            IFSObj.VC = [IFSObj.VC uint8(ones(size(IFSObj.VC,1),1)*255)];
        end;
        if size(IFSObj.VC,2) ~= 4
            error('Only [r,g,b] or [r,g,b,a] color format is supported.');
        end;
        
    % Handle faces colors:
    case {'fcol','fcolo','fcolor','fcolors','facecolor','facecolors'}
        if isempty(varargin{i+1}) || ~isnumeric(varargin{i+1})
            error('Face color value must be supplied as a matrix.');
        end;        
        IFSObj.FC = varargin{i+1};

        if ~isa(IFSObj.FC,'uint8')
            IFSObj.FC = uint8(255*IFSObj.FC);
        end;

        if size(IFSObj.FC,2) == 3
            IFSObj.FC = [IFSObj.FC uint8(ones(size(IFSObj.FC,1),1)*255)];
        end;
        if size(IFSObj.FC,2) ~= 4
            error('Only [r,g,b] or [r,g,b,a] color format is supported.');
        end;
                
    % Handling texture coordinates        
    case {'uv','tcoord','texturecoord','texturecoordinate','texturecoordinates'}
        if isempty(varargin{i+1}) || ~isnumeric(varargin{i+1})
            error('Texture-coordinates value must be supplied as a matrix.');
        end;        
        IFSObj.UV = double(varargin{i+1});
            
        if size(IFSObj.UV,2) > 3
            error('Maximum texture coordinates dimension is 3.');
        end;
    
    % Handling texture image
    case {'timg','timage','textureimage'}
        if isempty(varargin{i+1}) || ~isnumeric(varargin{i+1})
            error('Texture-image value must be supplied as a matrix.');
        end;        
        IFSObj.TImg = uint8(varargin{i+1});

    % Handling animation key-frames     
    case {'keyframe','keyframes'}
        if isempty(varargin{i+1}) || ~iscell(varargin{i+1})
            error('Animation key-frames value must be supplied as a cell array.');
        end;        
        A = varargin{i+1};
        IFSObj.A = cell(1,length(A));
        for a = 1:length(A),
            if size(A{a},2) > 3
                error('Maximum key-frame vertex dimension is 3.');
            end;
            
            if size(A{a},2) < 3
                A{a} = [A{a} zeros(size(A{a},1),3-size(A{a},2))];
            end;
                                    
            IFSObj.A{a} = single(A{a});
        end;
        clear a
        
    otherwise
        error([varargin{i} ' is not supported in IndexedFaceSet class.']);
    end;
end;
clear istart

% Checking IndexedFaceSet Validity:
if IFSObj.nv == 0 || IFSObj.nf == 0
    error('IndexedFaceSet geometry and connectivity are mandatory.');
end;

% Checking faces index validity:
mini = 1; maxi = IFSObj.nv;
for i=1:IFSObj.nf,
    mini = min(IFSObj.F{i});
    maxi = max(IFSObj.F{i});
end;
if mini < 1 || maxi > IFSObj.nv
    error('Face indices must be between 1 to number of vertices.');
end;
clear mini maxi

% Checking normals validity:
if ~isempty(IFSObj.VN) && size(IFSObj.VN,1) ~= IFSObj.nv
    error('Vertex normals and number of vertices must be the same.');
end;

if ~isempty(IFSObj.FN) && size(IFSObj.FN,1) ~= IFSObj.nf
    error('Face normals and number of faces must be the same.');
end;

% Checking color validity:
if ~isempty(IFSObj.VC) && size(IFSObj.VC,1) ~= IFSObj.nv
    error('Vertex colors and number of vertices must be the same.');
end;

if ~isempty(IFSObj.FC) && size(IFSObj.FC,1) ~= IFSObj.nf
    error('Face colors and number of faces must be the same.');
end;

% Checking texture-coordinates validity:
if ~isempty(IFSObj.UV) && size(IFSObj.UV,1) ~= IFSObj.nv
    error('Texture-coordinates and number of vertices must be the same.');
end;
    
% Checking key-frames animation:
if ~isempty(IFSObj.A)
    for a=1:length(IFSObj.A),
        if size(IFSObj.A{a},1) ~= IFSObj.nv
            error('Key-frames and number of vertices must be the same.');
		end;
    end;
    clear a
end;

IFSObj.nv = uint32(IFSObj.nv);
IFSObj.nf = uint32(IFSObj.nf);

% Constructing IndexedFaceSet class object
IFSObj = class(IFSObj,'IndexedFaceSet');
