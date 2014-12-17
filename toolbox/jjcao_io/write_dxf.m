function[]=write_dxf(fname,X,Y,Z) 
%
% Given a filename and a 3d mesh specified by X,Y, and Z 
% arrays it writes a DXF file with the surface specified 
% by the mesh, capable of being read by most CAD programs.
%
% writedxf.m by Greg Siegle

fullname=sprintf('%s.dxf',fname);
fid=fopen(fullname,'w');
fprintf(fid,'0\nSECTION\n2\nHEADER\n0\nENDSEC\n0\nSECTION\n2\nENTITIES\n0\n');
for a=1:size(X,1)-1
   for b=1:size(X,2)-1
      fprintf(fid,'3DFACE\n8\n0\n');
      %top left corner
      fprintf(fid,'10\n%.4f\n20\n%.4f\n30\n%.4f\n',X(a,b),Y(a,b),Z(a,b));
      %top right corner
      fprintf(fid,'11\n%.4f\n21\n%.4f\n31\n%.4f\n',X(a+1,b),Y(a+1,b),Z(a+1,b));      
      %bottom right corner
      fprintf(fid,'12\n%.4f\n22\n%.4f\n32\n%.4f\n',X(a+1,b+1),Y(a+1,b+1),Z(a+1,b+1));       
      %bottom left corner      
      fprintf(fid,'13\n%.4f\n23\n%.4f\n33\n%.4f\n',X(a,b+1),Y(a,b+1),Z(a,b+1));
      fprintf(fid,'0\n');
   end
end
fprintf(fid,'ENDSEC\n0\nEOF\n');
fclose(fid);
