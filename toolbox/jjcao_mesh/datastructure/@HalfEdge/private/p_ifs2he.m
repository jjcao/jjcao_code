function HEObj = p_ifs2he(IFSObj,HEObj)

nv = nov(IFSObj);
nf = nof(IFSObj);

% Setting number of edges:
HEObj.ne = sum(face_size(IFSObj));
HEObj.nb = 0;

% Setting orig, destination, next, previous edges and edge face:
HEObj.EO = zeros(1,HEObj.ne);
HEObj.ED = zeros(1,HEObj.ne);
HEObj.EF = zeros(1,HEObj.ne);
HEObj.EN = zeros(1,HEObj.ne);
HEObj.EP = zeros(1,HEObj.ne);

HEObj.VE = zeros(1,nv);
HEObj.FE = zeros(1,nf);

F = face_vertices(IFSObj);
if nof(IFSObj) == 1
    F = {F};
end

e = 1;
for f = 1:nf
    cf = F{f};
    ne = length(F{f});
    cf = [cf cf(1)];

    HEObj.EO(e:e+ne-1) = cf(1:end-1);
    HEObj.ED(e:e+ne-1) = cf(2:end);
    HEObj.EF(e:e+ne-1) = f*ones(1,ne);
    HEObj.EN(e:e+ne-1) = [e+1:e+ne-1 e];
    HEObj.EP(e:e+ne-1) = [e+ne-1 e:e+ne-2];

    HEObj.FE(f) = e;
    HEObj.VE(cf(1:end-1)) = e:e+ne-1;
    
    e = e+ne;
end;
clear F e f cf ne

% Setting edge twin and boundaries:
HEObj.ET = zeros(1,HEObj.ne);
EM = sparse(HEObj.EO,HEObj.ED,1:HEObj.ne,nv,nv);
[i,j] = find(EM ~= 0);
IJ = (j-1)*nv+i;
JI = (i-1)*nv+j;
HEObj.ET(EM(IJ)) = EM(JI);
clear EM i j IJ JI

% Adding boundary edges:
nte = find(HEObj.ET == 0);
HEObj.EO = [HEObj.EO HEObj.ED(nte)];
HEObj.ED = [HEObj.ED HEObj.EO(nte)];
HEObj.EF = [HEObj.EF zeros(1,length(nte))];
HEObj.ET = [HEObj.ET nte];
HEObj.ET(nte) = HEObj.ne+1:HEObj.ne+length(nte);
HEObj.ne = HEObj.ne + length(nte); 
HEObj.IBF = zeros(1,nf);
clear nte


% Tracing the boundaries:
ebl = find(HEObj.EF == 0);
if ~isempty(ebl) 
    cb = 1; nbe = ebl(1); pbe = nbe; cebl = nbe;    
    ebl = ebl(2:end);
    
    while ~isempty(ebl)
        nbe = ebl(find(HEObj.EO(ebl) == HEObj.ED(pbe)));
        if ~isempty(nbe)
            
            if length(nbe) > 1
                error(['Multiple bounday path (vertex #' num2str(HEObj.ED(pbe)) ')']);
            end;
            
            HEObj.EN(pbe) = nbe;
            HEObj.EP(nbe) = pbe;
            
            cebl = [cebl nbe];
            ebl = setdiff(ebl, nbe);            
            pbe = nbe;
            
        else
        
            if HEObj.EO(cebl(1)) ~= HEObj.ED(pbe)
                error('Boundary is not close');
            end;
            
            HEObj.EN(pbe) = cebl(1);
            HEObj.EP(cebl(1)) = pbe;
            HEObj.EF(cebl) = nf+cb;
            HEObj.BFS(cb) = length(cebl);
            HEObj.IBF(nf+cb) = 1;
            HEObj.FE(nf+cb) = cebl(1);

            cb = cb+1; nbe = ebl(1); pbe = nbe; cebl = nbe;
            ebl = setdiff(ebl,nbe);
            
        end;
    end;

    if HEObj.EO(cebl(1)) ~= HEObj.ED(pbe)
        error('Boundary is not close');
    end;
    
    HEObj.EN(pbe) = cebl(1);
    HEObj.EP(cebl(1)) = pbe;
    HEObj.EF(cebl) = nf+cb;
    HEObj.BFS(cb) = length(cebl);
    HEObj.IBF(nf+cb) = 1;
    HEObj.FE(nf+cb) = cebl(1);
    
    HEObj.nb = cb;

    clear ebl cb nbe pbe cebl
end;

% Setting the right data-type:
HEObj.ne = uint32(HEObj.ne);
HEObj.EO = uint32(HEObj.EO);
HEObj.ED = uint32(HEObj.ED);
HEObj.ET = uint32(HEObj.ET);
HEObj.EF = uint32(HEObj.EF);
HEObj.EN = uint32(HEObj.EN);
HEObj.EP = uint32(HEObj.EP);

HEObj.VE = uint32(HEObj.VE);
HEObj.FE = uint32(HEObj.FE);

HEObj.nb = uint32(HEObj.nb);
HEObj.BFS = uint32(HEObj.BFS);
HEObj.IBF = uint8(HEObj.IBF);
