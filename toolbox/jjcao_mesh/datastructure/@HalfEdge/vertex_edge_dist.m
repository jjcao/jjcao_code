function retval = vertex_edge_dist(HEObj,V)

ED = {};
for i=1:length(V),
   ED{i} = vedist(HEObj,V(i));
end

if length(ED) == 1
    retval = ED{1};
else
    retval = ED;
end




function VED = vedist(HEObj,v)

n = nov(HEObj);

D = inf*ones(1,n);
D(v) = 0;

Q = 1:n;
while ~isempty(Q)
    [d,u] = min(D(Q));
    
    Nu = vertex_neighbors(HEObj,Q(u));
    if isempty(Nu), return; end
    
    D(Nu) = min([D(Nu); (D(Q(u))+1)*ones(1,length(Nu))]);
    
    if u == 1
        Q = Q(2:end);
    elseif u == length(Q)
        Q = Q(1:end-1);
    else
        Q = [Q(1:u-1) Q(u+1:end)];
    end
end

VED = D;
