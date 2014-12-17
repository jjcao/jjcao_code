% Version: Nov. 2004
% Author : Chen Yanover
%
% sparse_cell: creates a sparse cell
%     S = sparse_cell(a) converts a sparse or full cell to sparse form by
%     squeezing out any empty elements.
%
%     S = sparse_cell(i,j,s,m,n) uses the rows of [i(:),j(:),s(:)] to generate
%     an m-by-n sparse cell.  The two integer index vectors, i and j, and the
%     entries cell array s, all have the same length.
%
%     There are several simplifications of this five argument call:
%     1. S = sparse_cell(i,j,s) uses m = max(i) and n = max(j).
%     2. S = sparse_cell(m,n) or sparse_cell([m,n]) abbreviates sparse_cell([],[],[],m,n).
%        This generates the ultimate sparse cell, an m-by-n all empty cell array.
%     3. S = sparse_cell creates an empty sparse_cell
% 
% Revised: March 2009; m,n not used any more, kept for backward compatibility 
function sc = sparse_cell(i,j,s,n,m)

% class members:
%   nelem   - number of non empty elements
%   indMat  - indices into the cell array
%   cells   - cell array
switch nargin,
  case 0,     % creates an empty sparse cell
    sc.nelem = 0;
    sc.indMat = sparse([]);
    sc.cells = {};
    sc = class(sc, 'sparse_cell');
  case 1,     % converts a sparse or full cell to sparse form	 
	 if isa(i, 'sparse_cell') sc=i;
    elseif iscell(i),
      [m, n] = size(i);
      nelem = 0;
      indMat = sparse(m, n);
      cells = {};
      for m_=1:m,
        for n_=1:n,
          if ~isempty(i{m_,n_}),
            nelem = nelem+1;
            indMat(m_,n_) = nelem;
            cells{end+1} = i{m_,n_};
          end
        end
      end
      sc.nelem = nelem;
      sc.indMat = indMat;
      sc.cells = cells;
      sc = class(sc, 'sparse_cell');
		elseif isnumeric(i) & prod(size(i))==2,
      sc = sparse_cell(i(1), i(2));
    else
      error('sparse_cell: wrong usage');
    end
  case {2,3,5},
    if nargin==2,
      m=i; n=j;
      i=[]; j=[]; s={};
    elseif nargin==3,
      m = max(i);
      n = max(j);
    elseif nargin==4, error('sparse_cell: wrong usage');
    end

    if ~iscell(s),
      error('sparse_cell: conversion to cell from double is not possible'); end

    i=i(:); j=j(:); s=s(:);
    if length(i)~=length(j) | length(j)~=length(s),
      error('sparse_cell: vectors must be the same lengths'); end

    sc.nelem = length(i);
    sc.indMat = sparse(i, j, 1:sc.nelem, m, n);
    sc.cells = s;
    sc = class(sc, 'sparse_cell');
  otherwise, error('sparse_cell: wrong usage');
end