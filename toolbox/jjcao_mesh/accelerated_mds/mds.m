% mds    Multidimensional scaling
%
% Usage:
%
%  X = mds(D,options)
%  [X,hist] = mds(D,options)
%
% Description: 
% 
%  Finds a configuration of points in Euclidean space with distances 
%  approximating the input matrix of distances D by solving a multidimensional 
%  scaling problem by SMACOF algorithm with vector extrapolation or multigrid acceleration. 
%
% Input:  
% 
%  D       - symmetric matrix of distances
%  options - structure of options created using set_options and containing
%            the following fields: 
%             options.X0      - initialization, matrix of size N*dim (default: randn(N,dim))
%             options.dim     - embedding dimension
%                               NOTE: either dim or x0 must be set 
%             options.iter    - number of iterations in SMACOF method/
%                               internal iterations in RRE method (default: 25)
%             options.cycles  - number of cycles in RRE method (default: 5)
%                               number of cycles in MG method (default: 3)
%             options.IND     - cell array with indices of points used at
%                               each grid, must be set in MG method
%             options.UPMTX   - cell array with interpolation matrices, must be set in MG method
%             options.DOWNMTX - cell array with decimation matrices, must be set in MG method
%             options.lambda  - parameter of modified stress in MG method (default: 1)
%             options.rtol    - relative stress change tolerance (default: 0.01)
%             options.atol    - absolute stress tolerance (default: 0)
%             options.method  - algorithm [smacof | rre | mg] (default: rre)
%             options.display - level of display [cycle | iter] 
%                               (default: iter)
%             options.xhistory- store solutions at each iteration [on | off] (default: off)
%
% Output:
%
%  X       - configuration of points in Euclidean space
%  hist    - history structure containing the following fields:
%             hist.s          - vector of stress values per iteration
%             hist.time       - vector of iteration duration in seconds
%             hist.X          - cell array of solutions at each iteration (optional)
%
% References:
%
%  [1] G. Rosman, A. M. Bronstein, M. M. Bronstein, A. Sidi, R. Kimmel,
%  "Fast multidimensional scaling using vector extrapolation", Techn. Report CIS-2008-01, 
%  Dept. of Computer Science, Technion, Israel, January 2008.
%
%  [2] G. Rosman, M. M. Bronstein, A. M. Bronstein, and R. Kimmel,
%  "Topologically constrained isometric embedding", Human Motion - Understanding, 
%  Modeling, Capture and Animation, Computational Imaging and Vision, Vol. 36, Springer, 2008.
%
%  [3] M. M. Bronstein, A. M. Bronstein, R. Kimmel, I. Yavneh, 
%  "Multigrid multidimensional scaling", Numerical Linear Algebra with Applications (NLAA), 
%  Special issue on multigrid methods, Vol. 13/2-3, pp. 149-171, March-April 2006.
%
%  [4] A. M. Bronstein, M. M. Bronstein, R. Kimmel, 
%  "Numerical geometry of nonrigid shapes", Springer, 2008.
%
% TOSCA = Toolbox for Surface Comparison and Analysis
% Web: http://tosca.cs.technion.ac.il
% Version: 1.1
%
% (C) Copyright Michael Bronstein, 2008
% All rights reserved.
%
% License:
%
% ANY ACADEMIC USE OF THIS CODE MUST CITE THE ABOVE REFERENCES. 
% ANY COMMERCIAL USE PROHIBITED. PLEASE CONTACT THE AUTHORS FOR 
% LICENSING TERMS. PROTECTED BY INTERNATIONAL INTELLECTUAL PROPERTY 
% LAWS. PATENTS PENDING.

function [X,hist] = mds(D,options)

% check input correctness
if size(D,1) ~= size(D,2),
    error('Matrix D must be square, exiting.')     
end

if strcmp(lower(options.method),'mg'),
    if ~isfield(options,'UPMTX') 
        error('No interpolation matrices provided, exiting') 
    end
    if ~isfield(options,'DOWNMTX') 
        error('No decimation matrices provided, exiting') 
    end
    if ~isfield(options,'IND') 
        error('No grid hierarchy provided, exiting') 
    end
    
    % set default
    if ~isfield(options,'lambda'),
        options.lambda = 1; 
    end
end

% initialization
if ~isfield(options,'X0')
   if isfield(options,'dim')
      options.X0 = randn(size(D,1),options.dim);
   else
      error('Unknown initialization/dimension, exiting. Set dim or X0.') 
   end
end

% set defauls
if ~isfield(options,'cycles'),
   options.cycles = 5; 
end

if ~isfield(options,'rtol'),
   options.rtol = 0.01; 
end

if ~isfield(options,'atol'),
   options.atol = 0; 
end

if ~isfield(options,'method'),
   options.method = 'rre'; 
end

if ~isfield(options,'iter'),
   switch lower(options.method), 
       case 'rre', 
           options.iter = 10; 
       case 'smacof',
           options.iter = 50; 
       case 'mg',
           options.iter = 3; 
   end
end

if ~isfield(options,'display'),
   options.display = 'iter'; 
end

if ~isfield(options,'xhistory'),
   options.xhistory = 'off'; 
end


% set flags
if nargout == 2,
    HISTORY = 1;
else
    HISTORY = 0;
end

switch lower(options.xhistory),
    case 'off',
        XHISTORY = 0;
    case 'on',
        XHISTORY = 1;
end


% print
if strcmp(lower(options.display),'iter') | strcmp(lower(options.display),'cycle')
    fprintf(1,'Solving MDS problem of size %dx%d\nEmbedding dimension: %d\nMethod: %s\n',...
            size(D,1),size(D,1),size(options.X0,2),options.method) 
end


% start optimization
switch lower(options.method),
   case 'rre'      % vector extrapolation using RRE
        if HISTORY,
            [X,hist] = smacof_rre(D,options.X0,options.cycles,options.iter,...
                                 options.display,options.xhistory,options.rtol,options.atol);
        else
            X = smacof_rre(D,options.X0,options.cycles,options.iter,...
                           options.display,options.xhistory,options.rtol,options.atol);
        end

    case 'mg'      % vector extrapolation using MG
        if HISTORY,
            [X,hist] = smacof_mg(D,options.X0,options.cycles,options.iter,...
                                 options.display,options.xhistory,options.rtol,options.atol,...
                                 options.IND,options.UPMTX,options.DOWNMTX,options.lambda);
        else
            X = smacof_mg(D,options.X0,options.cycles,options.iter,...
                           options.display,options.xhistory,options.rtol,options.atol,...
                                 options.IND,options.UPMTX,options.DOWNMTX,options.lambda);
        end

    case 'smacof'   % SMACOF
        if HISTORY,
            [X,hist] = smacof(D,options.X0,options.iter,...
                              options.display,options.xhistory,options.rtol,options.atol);
        else
            X = smacof(D,options.X0,options.iter,...
                       options.display,options.xhistory,options.rtol,options.atol);
        end
        
   otherwise
      error('Invalid method, exiting. Use method=smacof|rre|mg')
end

