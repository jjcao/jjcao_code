function [saliency] = log_spectral_saliency(verts, faces, options)
%
% singal scale saliency
%
% jjcao @ 2014
%
options.null = 0;
bDebug = getoptions(options, 'bDebug', 0);
bNormalize = getoptions(options, 'bNormalize', 0);
bSymmetrize = getoptions(options, 'bSymmetrize', 1);
n = getoptions(options, 'n', 9);
%% compute L
A = triangulation2adjacency(faces); % adjacency matrix 
ind = find(A>0);
[I, J] = ind2sub(size(A), ind);
dist2 = sum((verts(I,:) - verts(J,:)).^2, 2);
W = sparse(I,J, 1.0\(dist2+0.0000001) );
d = sum(W,2);

if bSymmetrize==1 && bNormalize==0
    L = diag(d) - W;
elseif bSymmetrize==1 && bNormalize==1    
    L = speye(M.nverts) - diag(d.^(-1/2)) * W * diag(d.^(-1/2));
elseif bSymmetrize==0 && bNormalize==1
    L = speye(M.nverts) - diag(d.^(-1)) * W;
else
    error('Does not work with bSymmetrize=0 and bNormalize=0');    
end
NW = diag(d.^(-1)) * W;

%% compute eigvector according to eigvalue in increasing order of magnitude, 
% [eigv,eigvalue] = eigs(L,M.nverts,'m'); %  matrix eigv whose columns are the corresponding eigenvectors
[eigv,eigvalue] = eig(full(L)); 
eigvalue = diag(eigvalue);
[Hf, index] = sort(eigvalue);
eigvalue = diag(Hf);
eigv = eigv(:, index);

if bDebug
    sum( sum(L - eigv*eigvalue*eigv'))
    sprintf('L is symmetric: %d, normalized: %d, vertices: %d, rank: %d', ...
        bSymmetrize, bNormalize, length(L), rank(full(L)))
end

%%
Hf(1) = 1; % Hf(2); % jjcao: pay attention!!
Lf = log(abs(Hf));
% local averaging filter
Af = filter(ones(n,1),n,Lf); % means filter Lf using Jn = ones(n,1)/n
Rf = abs(Lf - Af);
tmp = exp(Rf);
S = eigv*diag(tmp)*eigv'*W;
% S = eigv*diag(exp(Rf))*eigv'*NW;  % jjcao: pay attention!!

if bDebug
    h=figure; 
    outerpos = get(h,'OuterPosition');
    outerpos(4) = 800;
    set(h,'OuterPosition',outerpos);
    subplot(3,1,1);
    plot([1:length(Hf)], Hf, 'b-'); xlabel('Frequency index');ylabel('Laplaican');
    subplot(3,1,2);
    % plot([1:length(Hf)-1], Lf(2:end), 'g-'); xlabel('Frequency index');ylabel('Log Laplaican');
    plot([1:length(Hf)], Lf, 'g-'); xlabel('Frequency index');ylabel('Log Laplaican');
    subplot(3,1,3);
    plot([1:length(Hf)], Rf, 'g-'); xlabel('Frequency index');ylabel('spectral irregularity');
%     subplot(4,1,4);
%     plot([1:length(Hf)], tmp, 'g-'); xlabel('Frequency index');ylabel('exp spectral irregularity');
end

S = abs(S); % jjcao: pay attention
saliency = sum(S,2); 
% saliency = log(saliency);
% saliency = sum(S(:,2:end),2);
saliency = (saliency - min(saliency))/(max(saliency) - min(saliency));

if bDebug
    figure;
    trisurf(faces,verts(:,1),verts(:,2),verts(:,3), ...
    'FaceVertexCData', saliency, 'FaceColor','interp','edgecolor', 'none');
    axis off;colorbar;mouse3d; 
end  
