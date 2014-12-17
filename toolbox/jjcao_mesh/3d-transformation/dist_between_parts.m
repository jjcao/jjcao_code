function [mean_dist, max_dist] = dist_between_parts(v1, v2)
[mean_dist1, max_dist1] = distBetweenParts(v1, v2);
[mean_dist2, max_dist2] = distBetweenParts(v2, v1);
if mean_dist1<mean_dist2
    mean_dist = mean_dist1;
    max_dist = max_dist1;
else
    mean_dist = mean_dist2;
    max_dist = max_dist2;
end

function [mean_dist, max_dist] = distBetweenParts(v1, v2)
n1 = size(v1,1);
n2 = size(v2,1);
dists = zeros(n1,1);
for i = 1:n1
    minErr = 999999;
    pt1 = v1(i,:);
    for j = 1:n2
        pt2 = v2(j,:);
        tmp = sqrt( sum( (pt1-pt2).^2 ) );
        if tmp < minErr
            minErr = tmp;
        end
    end
    dists(i) = minErr;
end
mean_dist = mean(dists);
max_dist = max(dists);
