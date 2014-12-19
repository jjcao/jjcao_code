function [pts] = compute_vertex_normalization(pts, diameter)
% scale to unitBox and move to origin
% jjcao1231@gmail.com, Feb, 2012
    bbox = [min(pts(:,1)), min(pts(:,2)), min(pts(:,3)), max(pts(:,1)), max(pts(:,2)), max(pts(:,3))];
    c = (bbox(4:6)+bbox(1:3))*0.5;            
    pts = pts - repmat(c, size(pts,1), 1);
    s = diameter / max(bbox(4:6)-bbox(1:3));% make the bbox's diagnol = 1.6. %1.0, 1.6
    pts = pts*s;
end