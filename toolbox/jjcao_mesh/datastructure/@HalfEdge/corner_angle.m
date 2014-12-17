function ang = corner_angle(HEObj,edge1,edge2)

eps = 1e-6;

e1 = vertex_coords(HEObj,edge_orig(HEObj,edge1)) - vertex_coords(HEObj,edge_dest(HEObj,edge1));
e2 = vertex_coords(HEObj,edge_dest(HEObj,edge2)) - vertex_coords(HEObj,edge_orig(HEObj,edge2));
cp = cross(e2, e1); dp = dot(e2,e1);   

le1 = norm(e1,'fro');
le2 = norm(e2,'fro');
lcp = norm(cp,'fro');

if le1 == 0, le1 = eps; end;
if le2 == 0, le2 = eps; end;

ang = atan2(lcp/(le1*le2), dp/(le1*le2));

while ang < 0, ang = ang + 2*pi; end
