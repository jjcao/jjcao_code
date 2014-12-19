%
% This function returns a similarity matrix for two sets of
% signatures g1 and g2 based on the Euclidean distance
%
%   S = simmat_euclidean(g1, g2)
%
% -----------------------------------------------------------------
% (C) Richard (Hao) Zhang (2006)
%
function M = simmat_euclidean(g1, g2)

n1 = length(g1(:,1));
n2 = length(g2(:,1));

S = zeros(n1, n2);

for i = 1:n1
    for j = 1:n2
        S(i,j) = norm(g1(i,:) - g2(j,:));
    end
end

% Create structure
M.max_value = max(max(S));
M.value = S;
