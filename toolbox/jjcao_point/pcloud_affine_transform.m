% adapt a little from Andrea Tagliasacchi
function pts = pcloud_affine_transform( verts, T )
% T: a 4x4 affine transform matrix
% verts: input points
    pts = zeros( size(verts) );

    % a fat matrix containing a point in every row
    homogeneous_in = [ verts, ones(length(verts),1) ]';
    
    % apply the transform
    homogeneous_out = (T*homogeneous_in)';
    
    % revert homogeneous info
    pts(:,1) = homogeneous_out(:,1) ./ homogeneous_out(:,4);
    pts(:,2) = homogeneous_out(:,2) ./ homogeneous_out(:,4);
    pts(:,3) = homogeneous_out(:,3) ./ homogeneous_out(:,4);    
    