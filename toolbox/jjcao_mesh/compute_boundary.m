function boundary=compute_boundary(face, options)

% compute_boundary - compute the vertices on the boundary of a 3D mesh
%
%   boundary=compute_boundary(face);
%
%   changed to extract multi boundaries by jjcao 2013
%   Copyright (c) 2007 Gabriel Peyre

options.null = 0;
verb = getoptions(options, 'verb', 1);

if size(face,1)<size(face,2)
    face=face';
end

nvert=max(max(face));
nface=size(face,1);

A=sparse(nvert,nvert);
for i=1:nface
    if verb
        progressbar(i,nface);
    end
    f=face(i,:);
    A(f(1),f(2))=A(f(1),f(2))+1;
    A(f(1),f(3))=A(f(1),f(3))+1;
    A(f(3),f(2))=A(f(3),f(2))+1;
end
A=A+A';

boundary = {};
tag = zeros(nvert,1);
for i=1:nvert
    if tag(i), continue, end;
    u=find(A(i,:)==1);
    if ~isempty(u) && ~tag(u(1))
        border=[i u(1)];
        tag(i) = 1; tag(u(1)) = 1;
        s=border(2);
        j=2;
        while(j<=nvert)
            u=find(A(s,:)==1);
            if length(u)~=2
                warning('problem in border');
            end
            if u(1)==border(j-1)
                s=u(2);
            else
                s=u(1);
            end
            if s~=border(1)
                border=[border s];
                tag(s) = 1;
            else
                break;
            end
            j=j+1;
        end

        if j>nvert
            warning('problem in border');
        end
        boundary{end+1} = border;
    end
end




%%% OLD %%%
function v = compute_boundary_old(faces)

nvert = max(face(:));
ring = compute_vertex_ring( face );

% compute boundary
v = -1;
for i=1:nvert   % first find a starting vertex
    f = ring{i};
    if f(end)<0
        v = i;
        break;
    end
end
if v<0
    error('No boundary found.');
end
boundary = [v];
prev = -1;
while true
    f = ring{v};
    if f(end)>=0
        error('Problem in boundary');
    end
    if f(1)~=prev
        prev = v;
        v = f(1);
    else
        prev = v;
        v = f(end-1);
    end
    if ~isempty( find(boundary==v) )
        % we have reach the begining of the boundary
        if v~=boundary(1)
            warning('Begining and end of boundary doesn''t match.');
        else
            break;
        end
    end
    boundary = [boundary,v];
end