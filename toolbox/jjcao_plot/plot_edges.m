function h = plot_edges(edges, verts, sizee, color)

% plot_edges - plot a list of edges
%
%   edges: Nx2 matrix of indices
%   verts: Mx3 matrix of coordinates
%   color: a character, such as 'b';
%          or a vector [1 0 0]
%          or a vector (1 colour known per edge, not per vertex!!), such as color = (1:size(edges,1))'; 
%
%   changed from Darren Engwirda's code
%
%   Copyright (c) JJCAO 2012

if nargin<3
    sizee = 3;
end
if nargin<4
    color = 'b';
end

if size(color,1)<2
    h = patch('faces', edges, 'vertices', verts, 'edgecolor', color, 'LineWidth', sizee);axis equal;
else
    %% setup a new set of vertices/edges/colours with duplicate vertices so that each edge gets it's correct colour
    nnum = 0;
    pnew = zeros(2 * size(edges, 1), 3); %% new vertices
    enew = zeros(1 * size(edges, 1), 2); %% new edge indices
    cnew = zeros(2 * size(edges, 1), 1); %% new edge colours - via vertices
    for j = 1 : size(edges, 1)
        n1 = edges(j, 1); %% old edge indices
        n2 = edges(j, 2);
        enew(j, 1) = nnum + 1; %% new edge indicies into pnew
        enew(j, 2) = nnum + 2;
        pnew(nnum + 1, :) = verts(n1, :); %% create duplicate vertices
        pnew(nnum + 2, :) = verts(n2, :);
        cnew(nnum + 1) = color(j); %% map single edge colour onto both vertices
        cnew(nnum + 2) = color(j);
        nnum = nnum + 2;
    end

    %% Draw the set efficiently via patch
%     tic
    h = patch('faces', enew, 'vertices', pnew, 'facevertexcdata', cnew, ...
        'edgecolor', 'flat', 'facecolor', 'none', 'LineWidth', sizee);
%     plot(pnew(:,1), pnew(:,2), pnew(:,3), 'b.');
    axis equal;
%     toc
end

%% old 
% function h = plot_edges(edges, vertex, sizee, color)
% 
% % plot_edges - plot a list of edges
% %
% %   Where edges should be an Nx2 matrix of indices and points should be an Mx3 matrix of coordinates
% %
% %
% %   Copyright (c) 2004 Gabriel Peyr?
% %   Changed by JJCAO 2011, 2012
% 
% if nargin<3
%     sizee = 3;
% end
% if nargin<4
%     color = 'b';
% end
% 
% patch('faces', edges, 'vertices', points, 'edgecolor', 'b');
% 
% x = [ vertex(1,edges(1,:)); vertex(1,edges(2,:)) ];
% y = [ vertex(2,edges(1,:)); vertex(2,edges(2,:)) ];
% if size(vertex,1)==2
%     h = line(x,y, 'LineWidth', sizee, 'Color', color);
% elseif size(vertex,1)==3
%     z = [ vertex(3,edges(1,:)); vertex(3,edges(2,:)) ];
%     h = line(x,y,z, 'LineWidth', sizee, 'Color', color);
% else
%     error('Works only for 2D and 3D plots');    
% end