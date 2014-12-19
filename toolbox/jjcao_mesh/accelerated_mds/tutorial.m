% Demonstration of multidimensional scaling 
%
% References:
%
%  [1] A. M. Bronstein, M. M. Bronstein, R. Kimmel, 
%  "Numerical geometry of nonrigid shapes", Springer, 2008.
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

% load and display Swiss roll dataset
load swissroll

figure(1), trisurf(swiss.TRIV,swiss.X,swiss.Y,swiss.Z); axis image;
title('Swiss roll surface');
drawnow

% embed using SMACOF
disp('Embedding the intrinsic geometry of the Swiss roll swiss into R^3 using SMACOF (no acceleration)...')
options.X0 = [swiss.X,swiss.Y,swiss.Z];
options.method = 'smacof';
options.xhistory = 'on';
[X_smacof,hist_smacof] = mds(swiss.D,options);


pause

% show result at each iteration
for k = 1:length(hist_smacof.time),
    figure(2), trisurf(swiss.TRIV,hist_smacof.X{k}(:,1),hist_smacof.X{k}(:,2),hist_smacof.X{k}(:,3)); axis image;
    title(sprintf('Iteration %d, Stress=%g',k,hist_smacof.s(k)));
    drawnow 
    pause(0.25);
end


pause

% plot convergence
figure(3), semilogy(cumsum(hist_smacof.time),hist_smacof.s,'k');
xlabel('CPU time (sec)'); ylabel('stress');
title('Convergence speed');
