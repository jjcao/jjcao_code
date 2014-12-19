%
% This function returns a similarity matrix for two sets of
% signatures g1 and g2 based on the Chi-square distance
%
%   S = simmat_chisquare(g1, g2)
%
% -----------------------------------------------------------------
% (C) Richard (Hao) Zhang and Oliver van Kaick (2006, 2007)
%
function M = simmat_chisquare(g1, g2)

% Get length of contours
n1 = length(g1(:,1));
n2 = length(g2(:,1));

% Init similarity matrix
S = zeros(n1, n2);

% Compute similarity for all pars
for i = 1:n1
    for j = 1:n2
        % Get current descriptors
        h1 = g1(i,:);
        h1 = h1(:); % Make sure it's a vector
        h2 = g2(j,:);
        h2 = h2(:);
        % Apply Dirichlet smoothing
        h1 = h1 + 0.000000000001;
        h2 = h2 + 0.000000000001;
        % Compute Chi-square distance
        S(i,j) = sum(((h1 - h2).^2)./(h1 + h2));
    end
end

% Create structure
M.max_value = max(max(S));
M.value = S;
