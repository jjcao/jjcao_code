function [face_patch, patch_area,landmark]=supervertex_by_farthest_sampling(verts,faces,nbr_landmarks)
% 
% Input: verts
%        faces
%        nbr_landmark， number of supervertex
         
% Output: agd
%         landmark, the index of center vertex of each patch
%         patch_area, the area of each patch
%         face_patch, face_patch(i) = j means that the ith face is blong to patch j
% 
% Pipeline：
% 1. 计算采样点landmark
% 1.1 调用perform_farthest_point_sampling_mesh.m
% 2. 计算voronoi patch，得到每个点属于哪一个采样点Q
% 2.1 调用 computer_voronoi_mesh 该函数又调用
% 2.1.1 check_face_vertex.m 和 perform_fast_marching_mesh.m
% 3. 计算 face_patch，得到每个面属于那个采样点
% 4. 计算cellArea,得到每个采样点对应的cell的面积
%
% Copyright (c) 2013 Shuhua Li, Junjie Cao

DEBUG = false;

%% compute voronoi mesh
nverts = size(verts,1);
options.W = ones(nverts,1); % the speed function, for now constant speed to perform uniform remeshing
landmark = perform_farthest_point_sampling_mesh(verts,faces,[],nbr_landmarks,options);% perform the sampling of the surface
[Q,DQ, voronoi_edges] = compute_voronoi_mesh(verts,faces, landmark); %  i-th 顶点属于 Q(i)-th landmark

%% display voronoi patch
if DEBUG
    options.voronoi_edges = voronoi_edges;
    options.start_points = landmark;
    figure('name','voronoi');
    plot_fast_marching_mesh(verts,faces, Q(:,1), [], options); 
    view3d zoom;
end

%% 计算每个face属于哪一个cell
%i-th face属于 faceIndex(i)-th cell;若faceIndex(i)为0，则不属于任何一个cell
face_patch=zeros(size(faces,1),1);

%% 将有2-3个顶点属于b_i的face归入 cell_i
for i=1:nbr_landmarks           % cell_i
    vertIndex=find(Q(:,1)==i);  %属于cell_i 的顶点的编号，列向量
    % indexTemp:  nfaces*3 ，index1(s,t)=1,则顶点faces(s,t)在 cell_i中，否则为0
    indexTemp=zeros(size(faces));  
    for j=1:length(vertIndex)   
        indexTemp=indexTemp+(faces==vertIndex(j));
    end
    %index_3p： 3个顶点都在cell_i中的face编号,列向量
    index_3p=find(indexTemp(:,1).*indexTemp(:,2).*indexTemp(:,3)); 
    indexTemp(index_3p,:)=0;
    %index_2p: 2个顶点都在cell_i中的face编号,列向量
    index_2p=find(indexTemp(:,1).*indexTemp(:,2)+indexTemp(:,1).*indexTemp(:,3)+indexTemp(:,2).*indexTemp(:,3));
    face_patch(index_3p)=i;
    face_patch(index_2p)=i;
end

%% 3个顶点属于不同的landmark b的face,根据face的重心到3个顶点对应的采样点距离决定face属于哪个cell
index_1p=find(face_patch==0);    %3个顶点属于不同的landmark b的face编号
for i=1:size(index_1p)          %第index_1p(i)个face
    triIndex=faces(index_1p(i),:); %face的3个顶点编号 3*1
    triVerts=verts(triIndex,:);   %face的3个顶点坐标 3*3
    faceCenter=sum(triVerts)/3;   %face的重心坐标
    distTriVerts=zeros(3,1);              %face的3个顶点到重心的距离 3*1
    for j=1:3
    distTriVerts(j)=sum((triVerts(j,:)-faceCenter).^2);
    end
    dist=distTriVerts+DQ(triIndex,1);
    idx=find(dist==min(dist));    %face属于哪个顶点
    face_patch(index_1p(i))=Q(triIndex(idx)); %face属于哪个cell
end   

%% display cell
if DEBUG
    figure('Name','cell'); set(gcf,'color','white'); 
    options.face_vertex_color =face_patch;
    h = plot_mesh(verts, faces, options);
    colormap(jet(nbr_landmarks));     view3d rot; lighting none;    hold on;
%  display landmarks
    ms=25;
   for i=1:nbr_landmarks
       cellCenter= verts(landmark(i),:);
       h = plot3( cellCenter(1),cellCenter(2), cellCenter(3), 'r.');
       set(h, 'MarkerSize', ms);    
   end
end

%% compute area(b_i)
patch_area=zeros(nbr_landmarks,1);  %cell_i的面积为 patch_area(i)
for i=1:nbr_landmarks           % cell_i
    index2=find(face_patch==i);  %cell_i 中面的编号，列向量
    cellFaces=faces(index2,:);  %cell_i 中面
    A = cross(verts(cellFaces(:,2),:)- verts(cellFaces(:,1),:), verts(cellFaces(:,3),:)- verts(cellFaces(:,1),:)); %外积
    patch_area(i) = sum(0.5 * sqrt(A(:,1).^2+A(:,2).^2+A(:,3).^2));  %面积
end
