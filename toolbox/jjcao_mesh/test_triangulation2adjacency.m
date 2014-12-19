addpath ../jjcao_plot;
addpath ../jjcao_common;
clear all;close all;clc;

verts = [0 0 0; 1 0 0; 1 1 0; 0 1 0];
faces = [1 2 3; 3 1 4];
plot_mesh(verts, faces);

OA = [0 1 1 1; 
      1 0 1 0; 
      1 1 0 1;
      1 0 1 0];
A = triangulation2adjacency(faces);
assert( sum(sum(A-OA)) == 0);

OB = [0        1 1.4142 1; 
      1        0 1        0; 
      1.4142 1 0        1;
      1        0 1        0];
B = triangulation2adjacency(faces,verts);
assert( abs( sum(sum(B-OB)) ) < 0.001);
