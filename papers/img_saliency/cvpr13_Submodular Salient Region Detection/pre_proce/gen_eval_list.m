function [EvalSet aggloClID] = gen_eval_list( descSpCLR, locNet, spImg, ...
                AdjMat, NUM_EVAL_PNT) 
%
%	Generate Evaluation points per image
%   Source can be placed in only evaluation points.
%
% -------
% INPUTS:
%
%   descSpCLR           : Color descriptors of SPs ([M x d]) 
%                         where M: # ofSPs, d: feat dimension (ex. rgb: 3)
%   locNet              : intra-network of an image ([M x M]) 
%   spImg               : Superpixel image 
%   AdjMat              : Neighborhood connectivity of SPs in an
%                         image ([M x M])
%   NUM_EVAL_PNT        : number of evaluation points
%
% --------
% OUTPUTS:
%
%   EvalSet             : The SPs where submodularity will be
%                         checked (local SP indices). ([NUM_EVAL_PNT x 1])
%   aggloClID           : Agglomerative clustering result ([Mx1])
%
% -------------
%
% Gunhee Kim
%

% % --------------------------------------------------------------%
%     Parameters
% % --------------------------------------------------------------%

    EPS = 10^-6 ;
    BIGV = 10^3 ;
    SIM_THRES = 0.01 ;

    % In order to be a evaluation point, the cluster size should be
    % larger than this value. 
    CLUSTER_THRES = 2 ;

    nSP = size(AdjMat,1) ;

% % --------------------------------------------------------------%
%     Run Agglomerative clustering to get the evaluation set. 
% % --------------------------------------------------------------%

    % pairwise distance between SPs.
    % distSP = squareform(pdist(descSpCLR,'euclid')); 
    distSP = fn_dist_l2(descSpCLR') ; 
    % distSP(distSP==0) = EPS ;

    % Neighboring relations between SPs.
    AdjMat = AdjMat+AdjMat' ;
    % When applying agglomerative clustering
    % We only consider the neighbor pixles.
    distSPNgh = BIGV*ones(nSP) ;
    distSPNgh(AdjMat>0) = distSP(AdjMat>0) ;
    distSPNgh = distSPNgh(tril(ones(nSP,nSP),-1)>0) ;
    %distSPNgh(distSPNgh==0) = BIGV ;
    clear distSP ;

    % agglomerative clustering
    Z = linkage(distSPNgh','single'); 
    aggloClID = cluster(Z,'maxclust',NUM_EVAL_PNT);
    cntClID=histc(aggloClID,1:NUM_EVAL_PNT) ;
    sigClID = find(cntClID>=CLUSTER_THRES) ;
    NUM_EVAL_PNT = length(sigClID) ;


    % Get evaluation points.
    % first get the SPs along the image boundaries.
    bdSP = unique([spImg(1,:) spImg(end,:) spImg(:,1)' spImg(:,end)']) ;
    bdSP(bdSP==0) = [] ;
    bdLink = zeros(nSP,1) ; bdLink(bdSP) = 1 ;
    tmpMat = fn_normout([[locNet 2*bdLink] ; [ones(1,nSP) 0]]) ;
    for i=1:NUM_EVAL_PNT,
        cl_ind = find(aggloClID==sigClID(i)) ;
        % the most connected one in each cluster becomes the
        % evaluation point of the cluster. 
        %{
        pval = sum(fn_normout(locNet(cl_ind,cl_ind))) ;
        [max_val max_ind] = max(pval) ;
        EvalSet(i) = cl_ind(max_ind) ;
        %}
        ncl_ind = length(cl_ind) ;
        % expected number of visits before absorption.
        % absorption occurs at cluster/image boundaries
        tmpval = ones(1, ncl_ind) / (eye(ncl_ind) - tmpMat(cl_ind,cl_ind)) ;
        [max_val max_ind] = max(tmpval) ;
        EvalSet(i) = cl_ind(max_ind) ;
    end

    % ignore trivial positions further

    % run clustering from EvalSet
    %[cluster_id] = run_cluster_from_rankset(locNet, EvalSet, 1:NUM_EVAL_PNT) ;
    %sigClID = find(histc(cluster_id, 1:NUM_EVAL_PNT)>=CLUSTER_THRES) ;
    %EvalSet = EvalSet(sigClID) ;

    % sort evaluation points
    EvalSet = sort(EvalSet) ;
end

