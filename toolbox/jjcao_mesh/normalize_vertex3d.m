function [pts, newlen] = normalize_vertex3d(pts, diaglen)
%
% move to origin
% make length of bbox's diagonal = diaglen
%
% newlen: is new length of bbox's diagonal
% 
% jjcao 2014
%
bbox = [min(pts(:,1)), min(pts(:,2)), min(pts(:,3)), max(pts(:,1)), max(pts(:,2)), max(pts(:,3))];
c = (bbox(4:6)+bbox(1:3))*0.5;            
pts = pts - repmat(c, size(pts,1), 1);
s = diaglen / sqrt(sum((bbox(4:6)-bbox(1:3)).^2));
pts = pts*s;

%% test
bbox = [min(pts(:,1)), min(pts(:,2)), min(pts(:,3)), max(pts(:,1)), max(pts(:,2)), max(pts(:,3))];
newlen = sqrt(sum((bbox(4:6)-bbox(1:3)).^2));

