% listing = dir('*.cpp');
% for i=1:length(listing)
%     mex(listing(i).name)
% end

mex -largeArrayDims affinityic.cpp
mex -largeArrayDims cimgnbmap.cpp
mex -largeArrayDims mex_w_times_x_symmetric.cpp
mex -largeArrayDims sparsifyc.cpp
mex -largeArrayDims spmtimesd.cpp