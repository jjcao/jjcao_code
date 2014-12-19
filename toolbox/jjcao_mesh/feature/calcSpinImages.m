function spinImages = calcSpinImages(mesh, binSize, spinImageSize)

% spinImages = calcSpinImages(mesh, binSize, spinImageSize)
%
% This is a quick implementation of the spin images. A. Johnson and M. Hebert, 
% "Using Spin Images for Efficient Object Recognition in Cluttered 3D
% Scenes", IEEE TRANS. ON PATTERN ANALYSIS AND MACHINE INTELLIGENCE.
%
% calculates spin images of a mesh which are used for mesh registration and
% 3D object recognition. 
%
%
% Input Argument: mesh - data structure .. must have vertexNormals
%
% Output Argument: spinImages - data structure of spinIm
%
% Copyright : This code is written by Ajmal Saeed Mian {ajmal@csse.uwa.edu.au}
%              Computer Science, The University of Western Australia. The code
%              may be used, modified and distributed for research purposes with
%              acknowledgement of the author and inclusion this copyright information.
%
% Disclaimer : This code is provided as is without any warrantly.



for ii = 1 : length(mesh.vertices)
    
    spinIm = zeros(spinImageSize);
    Pt = mesh.vertices(ii,:);
    Normal = mesh.vertexNormals(ii,:);
    if isnan(Normal(1))
        continue;
    end
    for jj = 1 : length(mesh.vertices)
        thisPt = mesh.vertices(jj,:);
        thisNormal = mesh.vertexNormals(jj,:);
        if acos(Normal*thisNormal')>pi/3
            continue;
        end
        alpha = norm(cross((thisPt - Pt),Normal));
        beta = Normal(1)*(thisPt(1)-Pt(1)) + Normal(2)*(thisPt(2)-Pt(2)) + Normal(3)*(thisPt(3)-Pt(3));
        %beta = beta/sqrt(Normal(1)^2 + Normal(2)^2 + Normal(3)^2);
        % shifting the center up
         
        beta = beta - spinImageSize*binSize/2;
        if beta > 0 || beta < -spinImageSize*binSize
            continue;
        end        
        b = floor(-beta/binSize) + 1;
        
        a = floor(alpha/binSize) + 1;
        if a > spinImageSize
            continue;
        end
        spinIm(a,b) = spinIm(a,b) + 1;        
    end
    spinImages(ii).spinIm = spinIm;
end   