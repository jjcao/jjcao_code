%
% This function uses numerical labelling to show the correspondence
% K given between two contours Y1 and Y2:
%
%   viz_matching(Y1, Y2, K)
%
% Note that Y1 and Y2 may have different vertex counts. K is an (k x 2)
% correspondence vector. Specifically, the K(i,1)-th vertex of contour 
% 1 is matched with the K(i,2)-th vertex of contour 2.
%
% In our labelling, vertices of the first contour are labelled 
% consecutively in counterclockwise order. Then the label attached to
% a particular vertex v in contour 1 appears beside the vertex in 
% contour 2 that is corresponded by K. A vertex in contour 2 may be
% multiply labelled or unlabelled since the correspondence is not
% meant to be 1-1.
%
% -------------------------------------------------------------------
% (C) Richard (Hao) Zhang and Varun Jain and Oliver van Kaick (2006, 2007)
%
function viz_matching(Y1, Y2, K, varargin)

% Get globals
global contours;
global open_contour;
global feature_points;

% Set label shift
shift = 0.03;

% Order matching
K = sortrows(K);

% Check for additional parameters
display = 'labels';
if nargin > 3
    display = varargin{1};
end

positioning = 'centered';
if nargin > 4
    positioning = varargin{2};
end

% Get contour lengths
n1 = length(Y1(:,1));
n2 = length(Y2(:,1));

% Normalize shapes
Y1 = area_normalize(Y1);
Y2 = area_normalize(Y2);

% now to make sure that the two contours do not overlap, we translate
% them appropriately so that contour 1 is to the south-west of
% contour 2
x_max1 = max(Y1(:,1));
x_min1 = min(Y1(:,1));
y_max1 = max(Y1(:,2));
y_min1 = min(Y1(:,2));
x_max2 = max(Y2(:,1));
x_min2 = min(Y2(:,1));
y_max2 = max(Y2(:,2));
y_min2 = min(Y2(:,2));

if strcmp(positioning, 'displaced')
    Y1(:,1) = Y1(:,1) - x_min1;
    Y1(:,2) = Y1(:,2) - y_min1;
    Y2(:,1) = Y2(:,1) - x_min2 + x_max1 - x_min1;
    Y2(:,2) = Y2(:,2) - y_min2 + y_max1 - y_min1;
elseif strcmp(positioning, 'xdisplaced')
    Y1(:,1) = Y1(:,1) - x_min1;
    Y2(:,1) = Y2(:,1) - x_min2 + x_max1 - x_min1;
elseif strcmp(positioning, 'ydisplaced')
    Y1(:,2) = Y1(:,2) - y_min1;
    Y2(:,2) = Y2(:,2) - y_min2 + y_max1 - y_min1;
end

% draw the two contours now
if contours
    if open_contour
        set(plot(Y1(:,1), Y1(:,2), 'bo-'), 'LineWidth', [1.0]);
        hold on;
        set(plot(Y2(:,1), Y2(:,2), 'rx-'), 'LineWidth', [1.0]);
    else
        set(plot([Y1(:,1)' Y1(1,1)], [Y1(:,2)' Y1(1,2)], 'bo-'), 'LineWidth', [1.0]);
        hold on;
        set(plot([Y2(:,1)' Y2(1,1)], [Y2(:,2)' Y2(1,2)], 'rx-'), 'LineWidth', [1.0]);
    end
else
end

% label the first contour and assign label(s) to vertices in the second
% contour
if strcmp(display, 'labels')
    if isempty(feature_points)
        % initialize labeling for second contour
        for i=1:n2
            c2(i).label = '(';
        end

        % label the first contour and assign label(s) to vertices in the second
        % contour
        for i = 1:n1
            idx = find(K(:, 1) == i);
            if size(idx, 1) > 0
                if K(idx(1), 2) > 0
                    % write numerical labels for contour 1
                    set(text(Y1(i,1), Y1(i,2), num2str(i)), 'color', 'r');
                
                    % attach i to appropriate vertex in contour 2
                    c2(K(idx(1), 2)).label = strcat(c2(K(idx(1), 2)).label, num2str(i), ',');
                end
            end
        end

        for i=1:n2
            if length(c2(i).label) > 1
                c2(i).label(length(c2(i).label)) = ')';
                set(text(Y2(i,1), Y2(i,2), c2(i).label), 'color', 'b');
            end
        end        
    else
        % label the first contour and assign label(s) to vertices in the second
        % contour
        for i = 1:n1
            idx = find(K(:, 1) == i);
            if length(idx) > 0
                % write numerical labels for contour 1
                if length(find(feature_points == i)) > 0
                    current = find(feature_points == i);
                    h = text(Y1(idx(1),1)+shift, Y1(idx(1),2)+shift, num2str(current));
                    set(h, 'color', 'r');
                    set(h,'Interpreter','latex');
                    set(h, 'fontsize', 16);
                end
            end
        end

        % close the brackets in labels for contour 2 and display labels
        for i=1:n2
            idx = find(K(:, 2) == i);
            if length(idx) > 0
                j = K(idx(1), 1);
                if length(find(feature_points == j)) > 0
                    current = find(feature_points == j);
                    h = text(Y2(i,1)+shift, Y2(i,2)+shift, num2str(current(1)));
                    set(h, 'color', 'b');
                    set(h,'Interpreter','latex');
                    set(h, 'fontsize', 16);
                end
            end
        end
    end
end

% plot lines
if strcmp(display, 'lines')
    P = zeros(2, 2);
    for i = 1:n1
        idx = find(K(:, 1) == i);
        if length(idx) > 0
            if K(idx(1), 2) > 0
                P(1, :) = Y1(K(idx(1), 1), :);
                P(2, :) = Y2(K(idx(1), 2), :);
                plot(P(:, 1), P(:, 2), 'k-');
            end
        end
    end
end
