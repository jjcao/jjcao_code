function [LF,varargout] = vertex_local_frame(HEObj,V,varargin)
    
    n = length(V);
    LF = cell(n,1);
    
    if ((length(varargin) > 1) || (nargout > 2)),
        error(' ');
    end
    
    if (length(varargin) > 0),
        Vref = varargin{1};
    else
        Vref = zeros(n,1);
    end
    
    for i=1:n,
        v = V(i);
        vref = Vref(i);
        LF{i} = zeros(3,3);

        % Estimating the normal:
        if (vertex_normals_exist(HEObj)),
            ni = vertex_normals(HEObj,v);

        else
            vn = vertex_neighbors(HEObj,v);
            vnc = vertex_coords(HEObj,vn);
            vcp = cross(vnc,[vnc(2:end,:);vnc(1,:)],2);
            fna = sqrt(sum(vcp.^2,2));
            w = fna./sum(fna);
            
            if (face_normals_exist(HEObj))
                fn = edge_face(HEObj,vertex_edge_orig(HEObj,v));     
                fnn = face_normals(HEObj,fn);              
            else
                fnn = vcp ./ (fna*ones(1,3));
            end

            ni = sum(fnn.*(w*ones(1,3)));
            ni = ni./norm(ni,'fro');
        end
        
        % Finding the neighbor that is most close to the tagent plane:
        if (vref == 0)
            vn = vertex_neighbors(HEObj,v);
            vnc = vertex_coords(HEObj,vn);
            vncd = vnc - ones(length(vn),1)*vertex_coords(HEObj,v);
            vncdn = vncd ./ (sqrt(sum(vncd.^2,2))*ones(1,size(vncd,2)));
            vdp = dot(ones(length(vn),1)*ni,vncdn,2);

            [stam1,stam2] = min(vdp);
            vref = vn(stam2);
            clear stam1 stam2
        end
        
        uj = vertex_coords(HEObj,vref)-vertex_coords(HEObj,v);
        uj = uj./norm(uj,'fro');
            
        Vref(i) = vref;
        
        % Finding the third vector:
        wij = cross(ni,uj);
        wij = wij./norm(wij);

        % Fixining uj:
        uj = cross(wij,ni);
    
        LF{i}(:,1) = ni';
        LF{i}(:,2) = uj';
        LF{i}(:,3) = wij';
    end
    
    if (nargout > 1)
        varargout{1} = Vref;
    end
    