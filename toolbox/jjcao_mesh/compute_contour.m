function [v, contour]= compute_contour(vertices, faces, rings, boundary,sf,vid,eps)
% compute the parameter v coresponding to the argument u.
% the speed is almost as same as compute_contour_old.m 白写一场，引以为戒！
% jjcao, 2009

%% find start position and start vertice id
u = sf(vid);
u_boundary = sf(boundary);

cid = find(u_boundary - u>=0, 1, 'first');
if isempty(cid)
    warning('no contour for value = %f', u);
    v = []; contour = [];
    return;
end

pts = [];
if u == 0% || u == 1
    cid = boundary(cid);
    pts = vertices(cid,:);
else
    bid = cid-1;
    bid = boundary(bid);
    cid = boundary(cid);
    pts = midp(u,sf(bid),sf(cid),vertices(bid,:),vertices(cid,:));
end

%% compute all segments which's scalar == u
Ntriangles = size(faces,1);
segments = [];
for f=1:Ntriangles
    [p1,p2] = draw_eco_line(u,vertices(faces(f,:)',:), sf(faces(f,:)'));
    if ~isempty(p1)
        [id1, id2, pts] = add_pts(pts, p1', p2',eps);
        segments = [segments; [id1 id2]];
    end
end

if isempty(segments)
    warning('no contour for u = %f!', u);
    v = 0; contour = [];
    return;
end

%% connect segments to contour
edges = segments(end,:);
segments(end,:)=[];

while ~isempty(segments)
    ind = find(segments(:,1) == edges(1,1));
    if ~isempty(ind)
        edges = [segments(ind,2:-1:1); edges];
        segments(ind,:)=[];
    else    
        ind = find(segments(:,2) == edges(1,1));
        if ~isempty(ind)
            edges = [segments(ind,:); edges];
            segments(ind,:)=[];
        end
    end
    
    ind = find(segments(:,1) == edges(end,2));
    if ~isempty(ind)
        edges = [edges; segments(ind,:)];
        segments(ind,:)=[];
    else    
        ind = find(segments(:,2) == edges(end,2));
        if ~isempty(ind)
            edges = [edges; segments(ind,2:-1:1)];
            segments(ind,:)=[];
        end
    end
end

if edges(1,1) == 1
    sel = [edges(:,1); edges(end,2)];
else
    sel = [edges(end,2); edges(end:-1:1,1)];
end
contour = pts(sel,:);

%% compute v
vid_in_ordered = 0;
vertex = vertices(vid,:);
for i = 1:size(contour,1)
    if sqrt(sum((contour(i,:) - vertex).^2)) < eps
        vid_in_ordered = i;
        break;
    end
end

%  计算累加弦长
sel = [2:size(contour,1) 1];
D = cumsum( sqrt( sum( (contour(sel,:)-contour).^2, 2 ) ) );
D = [0; D(1:length(D)-1)];
d = D(end);
D = D./d;
v = D(vid_in_ordered);


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

function [id1 id2 pts] = add_pts(pts, p1, p2, eps)
%

id1 = []; id2 = [];
for i=1:size(pts,1)
   if sqrt( sum((pts(i,:)-p1).^2, 2) ) < eps
       id1 = i;
       break;
   end
end

for i=1:size(pts,1)
   if sqrt( sum((pts(i,:)-p2).^2, 2) ) < eps
       id2 = i;
       break;
   end
end

if isempty(id1)
    pts = [pts; p1];
    id1 = size(pts,1);
end

if isempty(id2)
    pts = [pts; p2];
    id2 = size(pts,1);
end
return;

function x = midp(v,y2,y1,x2,x1)
x = ((v - y1)/(y2 - y1))*(x2-x1) + x1;