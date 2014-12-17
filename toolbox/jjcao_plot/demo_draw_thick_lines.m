%
%    Demo code for 'mex_draw_thick_lines_on_img.cpp'.
%    
%    Draw multiple thick lines on the image.
%    You can freely set the thickness and color of each line.
%    
%    First of all, you need to compile and link the mex file
%    > mex mex_draw_thick_lines_on_img.cpp
%
%    Copyright, Gunhee Kim (gunhee@cs.cmu.edu)
%    Computer Science Department, Carnegie Mellon University, 
%    October 19 2009
%
%    Please email me if you find any bugs in the code. 
%

clear all; close all;

% if you've already compiled the mex-c file, skip the following line.
mex mex_draw_thick_lines_on_img.cpp -O

% load image
InputImg = imread('lena.jpg') ; 
%img = rgb2gray(img) ;

% number of lines to draw
NumLines = 10 ;
% max thickness
MaxThickness = 10 ;

[imr imc imch] = size(InputImg) ; 

% randomly select the coordinates of end points of lines
CoordPnt = rand(NumLines, 4);
CoordPnt(:,1:2) = round(imr*CoordPnt(:,1:2)) ;
CoordPnt(:,3:4) = round(imc*CoordPnt(:,3:4)) ;
CoordPnt(CoordPnt<1) = 1 ;

% randomly select the thickness of lines
Thickness = rand(NumLines,1) ;
Thickness = round(MaxThickness*Thickness) ;
Thickness = Thickness + MaxThickness ;

% randomly select the line colors
LineColor = rand(NumLines, imch);
LineColor = round(255*LineColor) ;

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% run the function 'mex_draw_thick_line_on_image.cpp'
%
% Inputs:
%      InputImg: Input image (Grayscale or Color)
%      CoordPnt: Coordinates of end points of lines [r1 r2 c1 c2] 
%          ( n x 4 double array, n: # of lines)
%      Thickness: Thickness of lines 
%          (The value is integer, but the type is double.)
%          (n x 1 double array, n: # of lines)
%      LineColor: Line colors (The values should be integers from 0 to 255)
%          (n x 4 double array, n: # of lines)
%          (Data type: double, ex. [255 255 255]: white )
%          (The channels of 'InputImg' and 'LineColor' should be consistent)
%
% Output: 
%      OutputImg: Output image (The same format with InputImg) 
%
OutputImg = mex_draw_thick_lines_on_img(InputImg, CoordPnt, Thickness, LineColor) ;
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% draw the result
figure;
subplot(1,2,1); imshow(InputImg);
subplot(1,2,2); imshow(OutputImg);


