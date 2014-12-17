function [potential] = get_potential_lineq_lu(cond_mat, lap_mat, ...
                        pos_pole, neg_pole, ilu_setup)
%
%	Compute Potential
%   Solution by incomplete LU decomposition ('ilu' function).
%   See the Page 21 of 'Random Walks and electric networks'
%       by  P.G. Doyle and J.L. Snell
%
% -------
% INPUTS:
%
%	cond_mat        : conductance network ([N x N]) 
%                     where N is the number of vertices
%	lap_mat         : laplace matrix ([N x N]) 
%                     if lap_mat is empty, compute it from 'cond_mat'.
%                     if lap_mat is given, the 'cond_mat' is not used here.
%   pos_pole        : nodes attached to sources
%   neg_pole        : nodes attached to ground
%   ilu_setup       : parameters for the 'ilu' function
%
% --------
% OUTPUTS:           
%   potential       : potential
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

% Parameters for 'ilu'
if nargin<5, 
    ilu_setup=[];
    ilu_setup.type = 'ilutp';
end


pos_pole = unique(pos_pole) ;
pos_pole = pos_pole(:)' ;
num_pos = length(pos_pole) ;
num_neg = length(neg_pole) ;

% Compute Laplace
if isempty(lap_mat),
    lap_mat = full(diag(sum(cond_mat,2)) - cond_mat) ;
end

% indices of unlabeled points
N = size(lap_mat,1) ;
uind = 1:N ;

uind([pos_pole neg_pole]) = [];

% Find RH
b = -sum(lap_mat(uind,pos_pole),2) ;
%b = -L(uind,pos_pole)*eye(num_pos) ;

% incomplete LU decomposition
[L U P] = ilu(sparse(lap_mat(uind,uind)), ilu_setup) ;
[p, foo] = find(P); p(p) = 1:length(p) ;

% Solve the linear system
pot = U\(L\(b(p,:))) ;

% Solve the linear system using 'linsolve'.
% pot = linsolve(full(U), linsolve(full(L), b(p,:), lt_opts), ut_opts) ;


% compute potential
potential = zeros(N,1) ; 
potential(uind) = pot ; 
potential(pos_pole) = 1 ;
potential(neg_pole) = 0 ;
