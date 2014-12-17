function [v iter] = get_potential_power(cond_net, pos_pole, neg_pole, v, ...
                                             maxiter, tol)
%
%	Compute Potential of each node. 
%   Use the power iteration algorithm
%   See the Page 21 of 'Random Walks and electric networks'
%       by  P.G. Doyle and J.L. Snell
%
% -------
% INPUTS:
%
%	cond_net        : conductance network
%   pos_pole        : nodes attached to sources
%   neg_pole        : nodes attached to ground
%   v               : initial potential (voltage)
%   maxiter         : max iteration
%   tol             : tolerance
%
% --------
% OUTPUTS:           
%   v               : potential
%
% -------------
%
% Written by Gunhee Kim @ CMU, 2011.
% All rights reserved.
%

% if neg_pole is not assigned, report an error.
if isempty(neg_pole)
    error('Negative pole should be given.');
end

N = size(cond_net,1) ;

% default parameters
if nargin<4,    v = zeros(N,1);   end
if nargin<5,    maxiter = 10000;        end
if nargin<6,    tol = 1e-5;             end

% if the matrix is not normalized.
if abs(sum(cond_net(:))-N)>tol*N,
    cond_net = fn_normout(cond_net);
end


% (CAUTION) 
% Do not need to inverse because my temperature is 
% the weighted sum of my neighbors.
% cond_net = cond_net' ;
delta = 1;
iter = 0;
v_pre = v ;
v(pos_pole) = 1 ;
v(neg_pole) = 0 ;
while (delta > tol && iter < maxiter)
    % power iteration
    v = cond_net * v ;

    % enforce boundary condition        
    v(pos_pole) = 1 ; v(neg_pole) = 0 ;

    delta = norm(v-v_pre,1) ;
    v_pre = v ;
    iter = iter + 1 ;
end
