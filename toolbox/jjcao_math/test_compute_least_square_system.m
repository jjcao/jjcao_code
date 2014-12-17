%test three lss solvers 
%Guangyuzhong 23/8/2013

%%
clear;clc;close all;
%MYTOOLBOXROOT='E:/jjcaolib/toolbox/';
MYTOOLBOXROOT='../';
addpath ([MYTOOLBOXROOT 'jjcao_common'])
addpath ([MYTOOLBOXROOT 'jjcao_math'])

%% get adj_mat
nn = 1000;
dn = 100;
fea = rand(nn,dn);
options.NeighborMode = 'KNN';
options.k = 7;
options.WeightMode = 'HeatKernel';
options.t = 1;
W = constructW(fea,options);

%% get L,b,constraint_id,constraint_value
L = diag(sum(W,2)) - W; 
b = zeros(size(L,1),1);

num_constraint_id = 10;
constraint_id =  unidrnd(nn,num_constraint_id,1)';
constraint_value = ones(num_constraint_id,1)';

%% solve Ax = b using solver 1 and 2
options.solver = 1; 
tic;
fid1 = compute_least_square_system(L, b, constraint_id, constraint_value, options);
toc;
options.solver = 2;
tic;
fid2 = compute_least_square_system(L, b, constraint_id, constraint_value, options);
toc;
options.solver = 3;
tic;
fid3 = compute_least_square_system(L, b, constraint_id, constraint_value, options);
toc;

%% compare the answer fid1 and fid2
figure('Name','difference btw solver1&solver2');plot(fid1-fid2);
figure('Name','difference btw solver1&solver3');plot(fid1-fid3);

diff1_2 = norm(fid1-fid2,'fro')
diff1_3 = norm(fid1-fid3,'fro')