function [V,E] = build_inner_graph(cont,ims,bSmoothCont)
if ~exist('bSmoothCont')	bSmoothCont = 0;	end
n_pt= size(cont,1);
X	= cont(:,1);
Y	= cont(:,2);
V	= cont;
msk		= ims ;
fg_mask = msk;
E	= build_graph_contour_C(X,Y,fg_mask,bSmoothCont);
E	= E';
figure;
disp_graph(V,E(:,1:2));	