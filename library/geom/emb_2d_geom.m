function [coord,num]=emb_2d_geom(iel,nx1,nx2,ny1,ny2,w1,s1,w2,h1,h2);
%forms the nodal coordinates and numbering for a 2-d
% slope of 8-node quadrilaterals. Nodes numbering in the y-direction,
% elements numbered in the x-direction.
zero=0;pt5=0.5;one=1.0;coord=0;num=0;
 nxe=nx1+nx2;
 nye=ny1+ny2;
 nt=nx1*ny1;
 nc=(3*nye+2)*nx1+2*ny1;
 facx=s1/ny1;
 facy=h1/ny1;
 facb=(w1+s1)/nx1;
 facs=zero;
 if (ny2~=0);facs=h2/ny2;end
 frh=zero;
 if (nx2~=0);frh=w2/nx2;end
 if (iel<=nt)
   iq=floor((iel-1)/nx1)+1;
   ip=iel-(iq-1)*nx1;
 else
   iq=floor((iel-nt-1)/nxe)+ny1+1;
   ip=iel-nt-(iq-ny1-1)*nxe;
 end
 if (ip<=nx1)
   num(1)=(ip-1)*(3*nye+2)+2*iq+1;
   num(2)=num(1)-1;
   num(3)=num(1)-2;
   num(4)=(ip-1)*(3*nye+2)+2*nye+iq+1;
   num(5)=ip*(3*nye+2)+2*iq-1;
   num(6)=num(5)+1;
   num(7)=num(5)+2;
   num(8)=num(4)+1;
   if(iq<=ny1)
     coord(1,1)=(ip-one)*(w1+iq*facx)/nx1;
     coord(3,1)=(ip-one)*(w1+(iq-1)*facx)/nx1;
     coord(5,1)=ip*(w1+(iq-1)*facx)/nx1;
     coord(7,1)=ip*(w1+iq*facx)/nx1;
     coord(1,2)=-iq*facy;
     coord(3,2)=-(iq-1)*facy;
     coord(5,2)=-(iq-1)*facy;
     coord(7,2)=-iq*facy;
   else
     coord(1,1)=(ip-one)*facb;
     coord(3,1)=(ip-one)*facb;
     coord(5,1)=ip*facb;
     coord(7,1)=ip*facb;
     coord(1,2)=-h1-(iq-ny1)*facs;
     coord(3,2)=-h1-(iq-ny1-1)*facs;
     coord(5,2)=-h1-(iq-ny1-1)*facs;
     coord(7,2)=-h1-(iq-ny1)*facs;
   end
 else
   num(1)=nc+(ip-nx1-1)*(3*ny2+2)+2*(iq-ny1)+1;
   num(2)=num(1)-1;
   num(3)=num(1)-2;
   num(4)=nc+(ip-nx1-1)*(3*ny2+2)+2*ny2+iq-ny1+1;
   num(5)=nc+(ip-nx1)*(3*ny2+2)+2*(iq-ny1)-1;
   num(6)=num(5)+1;
   num(7)=num(5)+2;
   num(8)=num(4)+1;
   coord(1,1)=w1+s1+(ip-nx1-1)*frh;
   coord(3,1)=coord(1,1);
   coord(5,1)=w1+s1+(ip-nx1)*frh;
   coord(7,1)=coord(5,1);
   coord(1,2)=-h1-(iq-ny1)*facs;
   coord(3,2)=-h1-(iq-ny1-1)*facs;
   coord(5,2)=-h1-(iq-ny1-1)*facs;
   coord(7,2)=-h1-(iq-ny1)*facs;
 end
 coord(2:2:6,:)=pt5*(coord(1:2:5,:)+coord(3:2:7,:));
 coord(8,:)=pt5*(coord(7,:)+coord(1,:));
num=num';
end

