
% compile the mex files, only needs to be done once
% mex mex_Pw_d.c
% mex mex_EMstep.c
% mex mex_logL.c


% load a test dataset
load 20News_w100.mat

% set some variables
Learn.Verbosity = 1;
Learn.Max_Iterations = 200; 
Learn.heldout = .1; % for tempered EM only, percentage of held out data
Learn.Min_Likelihood_Change = 1;
Learn.Folding_Iterations = 20; % for TEM only: number of fiolding
                               % in iterations

Learn.TEM = 1; %tempered or not tempered

% start pLSA
[Pw_z,Pz_d,Pd,Li,perp,beta] = pLSA(Xtrain,[],50,Learn);

% fold in test data
Pz_q = pLSA_EMfold(Xtest,Pw_z,[],25,beta);

Learn.TEM = 0;
[Pw_z2,Pz_d2,Pd2,Li] = pLSA(Xtrain,[],50,Learn);

Pq = sum(Xtest)./sum(Xtest(:));


% fold in test data
Pz_q2 = pLSA_EMfold(Xtest,Pw_z2,[],25,1);

% and output perplexity
p1 = pLSA_logL(Xtest,Pw_z,Pz_q,Pq);
p2 = pLSA_logL(Xtest,Pw_z2,Pz_q2,Pq);








