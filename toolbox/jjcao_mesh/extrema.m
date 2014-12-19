function [maximaIndex,minimaIndex] = extrema(facet,fun,k,percent)
%Seeking the local maxima and minima points of the function, the maxima is
%defined as: whose value is larger than the values at it's k-ring neighbor
%points
%%Input
%facet---the facet of the mesh
%fun---the function on the mesh
%k---the number of the rings
%%Output
%maximaIndex--the indexs of the local maxima points
%minmaIndex---the indexs of the local minima points

% change by jjcao 2013
%Hui Wang, Oct. 24, 2011, wanghui19841109@gmail.com

% W = graphAdjacencyMatrix(facet);
W = triangulation2adjacency(facet);

n = size(W,1);
if k > 1
  neighbourMatrix = spones(W^k) - speye(n);
elseif k == 1
  neighbourMatrix = W;
else
  error('k must be a negative integer!');
end

[row,col] = find(neighbourMatrix);
valPlus = fun(row) > fun(col);
flagMatrixPlus = sparse(row,col,valPlus);
valMinus = fun(row) < fun(col);
flagMatrixMinus = sparse(row,col,valMinus);

num = sum(neighbourMatrix,2);
numPlus = sum(flagMatrixPlus,2);
numMinus = sum(flagMatrixMinus,2);

maximaIndex = find(numPlus >= percent * num);
minimaIndex = find(numMinus >= percent * num);
