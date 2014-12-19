% test_read_mesh

% mesh to be validated
ov = [0 0 0; ...
    2 0 0; ...
    1 1.7321 0];
of = [1 2 3];

% read
filename = '../data/regular_triangle.off'
[v, f] = read_mesh(filename);

% vvalidate
assert(sum(sum(ov-v))==0);
assert(sum(sum(of-f))==0);

disp('test_read_mesh passed.');