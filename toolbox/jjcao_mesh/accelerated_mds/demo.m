% Demonstration of multidimensional scaling with reduced rank extrapolation (RRE) 
% and multigtid (MG) acceleration
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

pause;

% embed using SMACOF
disp('Embedding the intrinsic geometry of the Swiss roll swiss into R^3 using SMACOF (no acceleration)...')
options.X0 = [swiss.X,swiss.Y,swiss.Z];
options.method = 'smacof';
[X_smacof,hist_smacof] = mds(swiss.D,options);

pause;

% embed using SMACOF with RRE
disp('Embedding the intrinsic geometry of the Swiss roll swiss into R^3 using SMACOF with RRE acceleration...')
options.method = 'rre';
options.atol = hist_smacof.s(end);      % stop when reaching the same stress as SMACOF
[X_rre,hist_rre] = mds(swiss.D,options);

% load multigrid matrices
load mg_matrices

% embed using SMACOF with multigrid
disp('Embedding the intrinsic geometry of the Swiss roll swiss into R^3 using SMACOF with MG acceleration...')
options.method = 'mg';
options.DOWNMTX = DOWNMTX;
options.UPMTX = UPMTX;
options.IND = IND;
[X_mg,hist_mg] = mds(swiss.D,options);


% display flattened Swiss roll surface
figure(2), trisurf(swiss.TRIV,X_rre(:,1),X_rre(:,2),X_rre(:,3)); axis image;
title('Flattened Swiss roll surface');

pause;

% plot convergence
figure(3), semilogy(cumsum(hist_smacof.time),hist_smacof.s,'k:',...
                    cumsum(hist_rre.time),hist_rre.s,'r',...
                    cumsum(hist_mg.time),hist_mg.s,'b');
xlabel('CPU time (sec)'); ylabel('stress');
title('Convergence speed');
legend('SMACOF','RRE','MG');






