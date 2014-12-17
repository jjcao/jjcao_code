function [cluster_id,cluster_info,probMat] = compute_actual_clusters(A,rank_set,rank_set_val,options)
%
% todo, 用最后一个下降区域，或者用2，都要结合cluster的连通性
% Copyright (c) 2012 Junjie Cao
% DEBUG=1;
merge_option = getoptions(options,'merge_option',true);
nclusters = getoptions(options,'nclusters',6);
if merge_option==2    
    thresMerge = getoptions(options,'thresMerge',0.1);
end

switch merge_option
    case 0 % not merge
        [cluster_id probMat] = run_cluster_from_rankset(A, rank_set(1:nclusters), 1:nclusters);% run clustering from rank set.    
    case 1
        % method 1: 
        %     For example, cube_f300.off, when nclusters=19,
        %            diff(rank_set_val) =  50.2131   38.8015   37.4520   29.8319   29.6209    4.9748    4.8172    4.7276    4.5600    4.4437    4.3786    2.6479    2.4413  2.4391    2.4104    2.3785    2.3727    0.8354
        %            It is easy yo find 6 segments with large potentials & 1 segment followed.
        %     For example, cube_f300.off, when M.nsegments=9,
        %            diff(rank_set_val) =  50.2131   38.8015   37.4520   29.8319   29.6209    4.9748    4.8172    4.7276
        %            It is easy yo find 6 segments with large potentials & 1 segment followed.
        tmp = diff(rank_set_val);
        actual_rank_set = rank_set(tmp>mean(tmp(2:end)));
        actual_rank_set = rank_set(1:(length(actual_rank_set)+1));
        % actual_rank_set=rank_set(1:15);
        nclusters = length(actual_rank_set);
        [cluster_id probMat] = run_cluster_from_rankset(A, rank_set(1:nclusters), 1:nclusters);% run clustering from rank set.
    case 2
        % method 2
        lena = length(A);
        while true        
            [cluster_id probMat] = run_cluster_from_rankset(A, rank_set(1:nclusters), 1:nclusters);% run clustering from rank set.
            tmp = probMat;
            [C1,I] = max(tmp,[],2);
            linearInd = sub2ind(size(tmp), (1:lena)', I);
            tmp(linearInd)=0;
            C2 = max(tmp,[],2);
            if min(abs(C1-C2))<thresMerge
                nclusters = nclusters-1;
            else
                break;
            end
        end
end % end of switch

nclusters = max(cluster_id);
cluster_info = zeros(nclusters,1);
for j=1:nclusters
      cluster_info(j) = sum(cluster_id==j);
end

