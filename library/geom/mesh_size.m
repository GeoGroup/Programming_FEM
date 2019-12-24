function [ nels,nn ] = mesh_size( element,nod )
global nxe nye nze;
% This subroutine returns the number of elements (nels) and the number
% of nodes (nn) in a 2-d geometry-created mesh.
% nod 单元节点数
nels=0;nn=0;
if strcmp(element,'triangle')
   nels=nxe*nye*2;
   if (nod==3);nn=(nxe+1)*(nye+1);end;
   if (nod==6);nn=(2*nxe+1)*(2*nye+1);end;
   if (nod==10);nn=(3*nxe+1)*(3*nye+1);end;
   if (nod==15);nn=(4*nxe+1)*(4*nye+1);end;
elseif strcmp(element,'quadrilateral')
   nels=nxe*nye;
   if (nod==4);nn=(nxe+1)*(nye+1);end;
   if (nod==5);nn=(nxe+1)*(nye+1)+nxe*nye;end;
   if (nod==8);nn=(2*nxe+1)*(nye+1)+(nxe+1)*nye;end;
   if (nod==9);nn=(2*nxe+1)*(2*nye+1);end;
 elseif strcmp(element,'hexahedron')
   nels=nxe*nye*nze;
   if (nod==8);nn=(nxe+1)*(nye+1)*(nze+1);end;
   if (nod==14);nn=4*nxe*nye*nze+2*(nxe*nye+nye*nze+nze*nxe)+nxe+nye+nze+1;end;
   if (nod==20);nn=((2*nxe+1)*(nze+1)+(nxe+1)*nze)*(nye+1)+(nxe+1)*(nze+1)*nye;end;
end
end

