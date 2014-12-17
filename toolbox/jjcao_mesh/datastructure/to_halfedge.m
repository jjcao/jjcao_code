function [Mhe, Mifs]=to_halfedge(vertices, faces, varargin)
% varargin: vertices, faces [, attrib-name, attrib-val,...]
% Mhe: the halfedge class
% Mifs: the indexed face set
%
%   Copyright (c) 2008 JJCAO

for i=1:size(faces,1)
    F(i,1)={faces(i,:)};
end

ifscmd = 'indexedfaceset(vertices,F)';
if ~isempty(varargin)
    for i = 1:length(varargin)
        ifscmd = [ifscmd ',varargin{' num2str(i) '}'];
    end
    ifscmd = [ifscmd ')'];
end

Mifs = eval(ifscmd);
Mhe = halfedge(Mifs);