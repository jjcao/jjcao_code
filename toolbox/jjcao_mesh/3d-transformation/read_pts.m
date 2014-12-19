function verts = read_pts(filename)

fid=fopen(filename);
nverts = 0;
while ~feof(fid)
    tmp = fgetl(fid);
    nverts = nverts+1;
end
frewind(fid);
[A,cnt] = fscanf(fid,'%f, %f, %f', 3*nverts);
A = reshape(A, 3, length(A)/3);
verts = A';