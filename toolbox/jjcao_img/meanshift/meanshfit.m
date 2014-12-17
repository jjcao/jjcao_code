%GiG !
function [origDataPoints,dataPoints,clusterCentroids,pointsToClusters] = meanshift(numDimensions,useKNNToGetH)
        H  = 1; %if useKNNToGetH is false, H will be used else knn will be used to determine H
        threshold = 1e-3;
        numPoints = 300; 
        k = 100; % when we use knn this is the k that is used.

        numClasses = 0;
        pointsToCluster2 = [];
        centroids = [];

        function[dataPoints] = getDataPoints(numPoints,numDimensions)
                numClasses = randi([2,8],1); % generate a random number between 2 and 8
                dataPoints = [];
                curNumberOfPoints = 0;

                randomMeans = [];
                for i = 1:numClasses
                        randomMean = randi([0,100], 1) * rand(1);
                        randomStd = randi([1,2]) * rand(1);
                        curGaussianModelPoints = randomMean + randomStd .* randn(numPoints,numDimensions);
                        dataPoints = [dataPoints;curGaussianModelPoints];
                        randomMeans = [randomMeans,randomMean];
                        pointsToCluster2(curNumberOfPoints +1 : (curNumberOfPoints+numPoints)) = i;
                        curNumberOfPoints = curNumberOfPoints + numPoints;
                        centroids = randomMeans;
                end
        end

        function d = sqdist(a,b)
                %taken from demo code of class
                aa = sum(a.*a,1); bb = sum(b.*b,1); ab = a'*b;
                d = abs(repmat(aa',[1 size(bb,2)]) + repmat(bb,[size(aa,2) 1]) - 2*ab);
        end

        function bandwidth = getBandWith(curDataPoint,origDataPoints,euclideanDist,useKNNToGetH)
                if (useKNNToGetH == false)
                        bandwidth = H;
                else
                        [sortedEuclideanDist,indexInOrig] = sort(euclideanDist);
                        bandwidth = sortedEuclideanDist(k);
                end
        end

        function [clusterCentroids,pointsToCluster] = getClusters(origDataPoints,finalDataPoints)
                [numSamples,numFeatures] = size(finalDataPoints);
                clusterCentroids  = [];
                pointsToCluster = [];
                numClusters = 0;

                for i = 1:numSamples
                        closestCentroid = 0 ;
                        curDataPoint = finalDataPoints(i,:);

                        for j = 1:numClusters
                                distToCentroid = sqrt(sqdist(curDataPoint',clusterCentroids(j,:)'));
                                %distToCentroid = sqdist(curDataPoint',clusterCentroids');
                                %if (distToCentroid <  8 * H) 
                if (distToCentroid <  4 * H) 
                                        closestCentroid = j;
                                        break;
                                end
                        end
                        if (closestCentroid > 0)
                                pointsToCluster(i,:) = closestCentroid;
                                clusterCentroids(closestCentroid,:) = 0.5 * (curDataPoint + clusterCentroids(closestCentroid,:));
                        else
                                numClusters = numClusters + 1 ;
                                clusterCentroids(numClusters,:) = finalDataPoints(i,:);
                                pointsToCluster(i,:) = numClusters;
                        end
                end
        end

        function plotPoints(clusterCentroids,pointsToCluster,origDataPoints)
                [numSamples,numFeatures] = size(origDataPoints);
                allColorsInMatlab = 'bgrcmyk'; %ignoring white as it is the default bg for the plot
                sizeOfColors = size(allColorsInMatlab,2);

                if (numFeatures == 2)
                        % plot the original Data Points
                        h = figure(1);
                        hold on;
                        for i=1:numClasses
                                colourIndex = mod(i,sizeOfColors) + 1;
                                allElemsInThisCluster = find(pointsToCluster2 == i);
                                allOrigPointsInThisCluster = origDataPoints(allElemsInThisCluster,:);
                                plot(allOrigPointsInThisCluster(:,1), allOrigPointsInThisCluster(:,2), [allColorsInMatlab(colourIndex) '.']);
                        end
                        plot(centroids,centroids,'s');
                        hold off;

                        h = figure(2);
                        hold on;
                        % plot the original Data Points
                        for i = 1:size(clusterCentroids,1) 
                                colourIndex = mod(i,sizeOfColors) + 1 ;
                                allElemsInThisCluster = find(pointsToCluster == i);
                                allOrigPointsInThisCluster = origDataPoints(allElemsInThisCluster,:);
                                plot(allOrigPointsInThisCluster(:,1), allOrigPointsInThisCluster(:,2), [allColorsInMatlab(colourIndex) '.']);
                        end
                        hold off;
        end
        end

        function [origDataPoints,dataPoints] = doMeanShift(dataPoints,useKNNToGetH)
                [numSamples,numFeatures] = size(dataPoints);

                origDataPoints = dataPoints;

                for i = 1:numSamples
                        diffBetweenIterations = 10;

                        while (diffBetweenIterations > threshold)
                                curDataPoint = dataPoints(i,:);
                                euclideanDist = sqdist(curDataPoint',origDataPoints');
                                bandwidth = getBandWith(origDataPoints(i,:),origDataPoints,euclideanDist,useKNNToGetH);
                                kernelDist = exp(-euclideanDist ./ (bandwidth^2));
                                numerator = kernelDist * origDataPoints;
                                denominator = sum(kernelDist);
                                newDataPoint = numerator/denominator;
                                dataPoints(i,:) = newDataPoint;
                                diffBetweenIterations = abs(curDataPoint - newDataPoint);
                        end
                end
                
                [clusterCentroids,pointsToClusters] = getClusters(origDataPoints,dataPoints);
                plotPoints(clusterCentroids,pointsToClusters,origDataPoints);
        end
        

        dataPoints = getDataPoints(numPoints,numDimensions);
        [origDataPoints,dataPoints] = doMeanShift(dataPoints,useKNNToGetH);
end
