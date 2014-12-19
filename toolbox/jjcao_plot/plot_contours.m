function [handle]=plot_contours(triang_faces,vertices,sf, options)
% [handle]=plot_contours(triang_faces,vertices,sf,Nlines);
% Draws isolines (also called contour lines) of a scalar function linearly
% interpolated on a triangulated surface by its sf at vertices.
% NOTATION:
% 'Ntriangles' denotes the number of triangles, which is the number
% of rows in the (input) array 'triang_faces'. 'Ngrid_pts' denotes the
% number of vertices which is the number of rows in the (input)
% array 'vertices'. 'DIM' is the dimension of the physical space,
% also the number of columns in 'vertices' (DIM==2 or DIM==3).
% PARAMETERS:
% -- triang_faces (Ntriangles x 3)
%     triang_faces(i,j) is the index of
%     the vertex no. 'j' (j=1,2,3) of the triangle no. 'i' into the
%     array 'vertices'.
% -- vertices (Ngrid_pts x DIM)
%     vertices(k,:) are the coordinates of the vertex no. 'k'.
% -- sf (Ngrid_pts x 1)
%     Value of the interpolated function, at each vertex.
% -- options
%    -- Nlines (1x1)
%       Number of levels of isolines to draw.
%    -- values
%       array for value of isolines to draw. If values are assigned, then
%    Nlines is just ingored.
%    -- lineWidth
% Returns handle to the graphics object.
%
% Copyright (c) 2009 jjcao, adapt from Dima Sorkin

lineWidth = getoptions(options, 'lineWidth', 2);
S.EdgeColor = getoptions(options, 'edgecolor', 'flat'); % 'black';
S.LineWidth = lineWidth;
S.CDataMapping = 'scaled';
S.faces = [];
S.vertices = [];

S.FaceVertexCData = [];

MAXA = max(sf(:));
MINA = min(sf(:));

values = getoptions(options, 'values', []);
if ~isempty(values)
    LV = values;
else
    Nlines = getoptions(options, 'Nlines', 10);
    LV = linspace(MINA,MAXA,Nlines);
end

v_ix = 1;

Ntriangles = size(triang_faces,1);
[Ngrid_pts,DIM] = size(vertices);

if Ngrid_pts ~= length(sf)
   error('Ngrid_pts ~= length(sf)');
end

for f=1:1:Ntriangles
    for l_ix= 1:1:length(LV)
        [p1,p2] = draw_eco_line(LV(l_ix),vertices(triang_faces(f,:)',:),...
                                sf(triang_faces(f,:)'));
        if ~isempty(p1)
           S.FaceVertexCData = [S.FaceVertexCData;LV(l_ix);LV(l_ix)];
           S.vertices = [S.vertices ; p1' ; p2'];
           S.faces = [S.faces ;[v_ix,v_ix + 1]];
           v_ix = v_ix +2 ;
        end
    end
end

handle = patch_checked(S);

grid on;
xlabel('x');
ylabel('y');
if DIM == 3
   zlabel('z');
end
colorbar;

function [posA,posB]=draw_eco_line(val,vertices,sf)
%

if (val < min(sf))  || (val > max(sf))
    posA = [] ; posB = [];
    return;
end

pos = [];

if (val >= min(sf(1:2)))  && (val <=  max(sf(1:2)))
    if (val ~= sf(2))
        pos = [pos;midp(val,sf(1),sf(2),vertices(1,:),vertices(2,:))];
    end
end

if (val >= min(sf(2:3)))  && (val <=  max(sf(2:3)))
    if (val ~= sf(3))
        pos = [pos;midp(val,sf(2),sf(3),vertices(2,:),vertices(3,:))];
    end
end

if (val >= min(sf([1,3])))  && (val <=  max(sf([1,3])))
    if(val ~= sf(1))
        pos = [pos;midp(val,sf(1),sf(3),vertices(1,:),vertices(3,:))];
    end
end

Npos = size(pos,1);
if Npos ~= 2
    tmp = find(sf==val);
    if (length(tmp)==2)
        pos = vertices(tmp,:);
    else
        posA=[];  posB=[];
        return;
    end    
end

posA = pos(1,:)';
posB = pos(2,:)';


function x = midp(v,y2,y1,x2,x1)
x = ((v - y1)/(y2 - y1))*(x2-x1) + x1;
