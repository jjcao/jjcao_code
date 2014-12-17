% extrap   Reduced rank exrapolation (RRE)
%
% Usage:
%
%  called from smacof_rre
%
% TOSCA = Toolbox for Surface Comparison and Analysis
% Web: http://tosca.cs.technion.ac.il
% Version: 1.0
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
function x_ = extrap(XX)
K = size(XX,2)-1;
dX  = XX(:,2:end) - XX(:,1:end-1);
[Q,R] = mgs(dX);
y = R'\ones(K,1);
gamma = R\y;
gamma = gamma(:)/sum(gamma);

x_ = sum(repmat(gamma',[size(XX,1) 1]).*XX(:,1:K),2);
return


% mr    Modified Gram-Schmidt. This is a more stable way to compute a
%       QR factorization
function [Q, R] = mgs(A);
[m,n] = size(A);
% we assume that m>= n.
V = A;
R = zeros(n,n);
for i=1:n,
   R(i,i) = norm(V(:,i));
   V(:,i) = V(:,i)/R(i,i);
   if (i < n)
      for j = i+1:n
         R(i,j) = V(:,i)' * V(:,j);
         V(:,j) = V(:,j) - R(i,j) * V(:,i);
      end
   end
end
Q = V;

