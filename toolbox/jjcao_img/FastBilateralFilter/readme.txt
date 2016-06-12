This is the MATLAB implementation of the fast O(1) 
bilateral filter described in the following papers:

[1] K.N. Chaudhury, D. Sage, and M. Unser, "Fast O(1) bilateral filtering using trigonometric range kernels," IEEE Trans. Image 
Processing, vol. 20, no. 11, 2011.

[2] K.N. Chaudhury, "Acceleration of the shiftable O(1) algorithm for
bilateral filtering and non-local means,"  IEEE Transactions on Image Proc., 
vol. 22, no. 4, 2013.

Also included is the MATLAB implementation of joint bilateral filtering 
of multiband images.

Author   : Kunal N. Chaudhury (kunal@ee.iisc.ernet.in)
Date       : July 2015

To run the software use:
=======================

[ f , param ]  =  shiftableBF(f0, sigma1, sigma2, w, tol)

INPUT
=====

f0              :  input grayscale image
f                : bilateral filter output
sigmas       : width of spatial Gaussian
sigmar       : width of range Gaussian
[-w, w]^2   : domain of spatial Gaussian
tol              : truncation error (kernel)


OUTPUT
======

f  : filtered output image
param  : list of parameters


 
