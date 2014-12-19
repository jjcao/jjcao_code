% Added by Hoi Wong Mar/10/2009
function outputIndex = end(A, K, N)

% I don't think that's necessary, as an internal matlab function, it's probably never invoked
if K>N, 
	error('The index position cannot be greater than # of indices involved');
end

switch N, 
 case 1, outputIndex=prod(size(A));
 case 2, 
	if     K==1, outputIndex=size(A,1);
	elseif K==2, outputIndex=size(A,2);
	end
 otherwise, error('ND sparce cell is not supported'); 
end

