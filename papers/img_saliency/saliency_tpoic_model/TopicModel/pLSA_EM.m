% function [Pw_z,Pd_z,Pz,Li] = pLSA_EM(X,K,Par)
%
% Probabilistic Latent semantic alnalysis (pLSA)
%
% Notation:
% X ... (m x nd) term-document matrix (observed data)
%       X(i,j) stores number of occurrences of word i in document j
%
% m  ... number of words (vocabulary size)
% nd ... number of documents
% K  ... number of topics
%
% Li   ... likelihood for each iteration
% Pz   ... P(z)
% Pd_z ... P(d|z) 
% Pw_z ... P(w|z) corresponds to beta parameter in LDA
%
% Pz_wd ... P(z|w,d) posterior on z
%
% 
% References: 
% [1] Thomas Hofmann: Probabilistic Latent Semantic Analysis, 
% Proc. of the 15th Conf. on Uncertainty in Artificial Intelligence (UAI'99) 
% [2] Thomas Hofmann: Unsupervised Learning by Probabilistic Latent Semantic
% Analysis, Machine Learning Journal, 42(1), 2001, pp.177.196 
%
% Josef Sivic
% josef@robots.ox.ac.uk
% 30/7/2004

function [Pw_z,Pd_z,Pz,Li] = pLSA_EM(X,K,Par)

if nargin<3
   Par.maxit  = 100;
   Par.Leps   = 1;   
   Par.doplot = 0;
end;   

if Par.doplot
   ff(1)=figure(1); clf;
   %set(ff(1),'Position',[1137         772         452         344]);
   title('Log-likelihood');
   xlabel('Iteration'); ylabel('Log-likelihood');
   %figure(2); clf; 
end;


m  = size(X,1); % vocabulary size
nd = size(X,2); % # of documents

% initialize Pz, Pd_z,Pw_z
[Pz,Pd_z,Pw_z] = pLSA_init(m,nd,K);

% allocate memory for the posterior
Pz_dw = zeros(m,nd,K); 

Li    = [];
maxit = Par.maxit;

% EM algorithm
for it = 1:maxit   
   fprintf('Iteration %d ',it);
   
   % E-step
   Pz_dw = pLSA_Estep(Pw_z,Pd_z,Pz);
   
   % M-step
   [Pw_z,Pd_z,Pz] = pLSA_Mstep(X,Pz_dw);
   if Par.doplot>=2
     Pw_z
   end;  
   
   % Evaluate data log-likelihood
   Li(it) = pLSA_logL(X,Pw_z,Pz,Pd_z);   
        
   % plot loglikelihood
   if Par.doplot>=3
      figure(ff(1));
      plot(Li,'b.-');
   end;
      
   
   % convergence?
   dLi = 0;
   if it > 1
     dLi    = Li(it) - Li(it-1);
     if dLi < Par.Leps, break; end;   
   end;
   fprintf('dLi=%f \n',dLi);
end;



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% initialize conditional probabilities for EM 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [Pz,Pd_z,Pw_z] = pLSA_init(m,nd,K)
% m  ... number of words (vocabulary size)
% nd ... number of documents
% K  ... number of topics
%
% Pz   ... P(z)
% Pd_z ... P(d|z)
% Pw_z ... P(w|z)

Pz   = ones(K,1)/K; % uniform prior on topics

% random assignment
Pd_z = rand(nd,K);   % word probabilities conditioned on topic
C    = 1./sum(Pd_z,1);  % normalize to sum to 1
Pd_z = Pd_z * diag(C);

% random assignment
Pw_z = rand(m,K);
C    = 1./sum(Pw_z,1);    % normalize to sum to 1
Pw_z = Pw_z * diag(C);

return;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% (1) E step compute posterior on z,  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function Pz_dw = pLSA_Estep(Pw_z,Pd_z,Pz)
   K = size(Pw_z,2);

   for k = 1:K
      Pz_dw(:,:,k) = Pw_z(:,k) * Pd_z(:,k)' * Pz(k);
   end;   
   C = sum(Pz_dw,3);

   % normalize posterior
   for k = 1:K
      Pz_dw(:,:,k) = Pz_dw(:,:,k) .* (1./C);
   end;   
return;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  (2) M step, maximazize log-likelihood
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [Pw_z,Pd_z,Pz] = pLSA_Mstep(X,Pz_dw)
   K = size(Pz_dw,3);
   
   for k = 1:K
      Pw_z(:,k) = sum(X .* Pz_dw(:,:,k),2);
   end;
   
   for k = 1:K
      Pd_z(:,k) = sum(X.* Pz_dw(:,:,k),1)';
   end;
   
   Pz = sum(Pd_z,1);
   
   % normalize to sum to 1
   C = sum(Pw_z,1);
   Pw_z = Pw_z * diag(1./C);
   
   C = sum(Pd_z,1);
   Pd_z = Pd_z * diag(1./C);
   
   C = sum(Pz,2);
   Pz = Pz .* 1./C;
   Pz;
return;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% data log-likelihood
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function L = pLSA_logL(X,Pw_z,Pz,Pd_z)
   L = sum(sum(X .* log(Pw_z * diag(Pz) * Pd_z' + eps)));
return;



return;