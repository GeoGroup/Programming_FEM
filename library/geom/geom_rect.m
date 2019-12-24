function [ num,coord ] = geom_rect(iel,element,x_coords,y_coords,nod,dir )
% This subroutine forms the coordinates and connectivity for a
% rectangular mesh of rt. angled triangular elements (3, 6, 10 or 15-node)
% or quadrilateral elements (4, 8 or 9-node) counting in the
% x- or y-dir. 
jel=0;ip=0;iq=0;
pt5=0.5;two=2.0;d3=3.0; 
num=zeros(nod,1);
nxe=size(x_coords,1)-1;
coord=zeros(nod,2);
if strcmp(element,'triangle')
    nye=(size(y_coords,1)-1)*2;
    if strcmp(dir,'x') || strcmp(dir,'r')
        jel=2*nxe*floor((iel-1)/(2*nxe));
        ip=floor((iel-jel+1)/2);
        iq=2*(floor((iel-1)/(2*nxe))+1)-1+floor((floor(iel/2)*2)/iel);
    else
        jel=floor((iel-1)/nye);
        ip=jel+1;
        iq=iel-nye*jel;
    end
    switch nod
        case 3
            if mod(iq,2)~=0
                if strcmp(dir,'x') || strcmp(dir,'r')
                    num(1)=floor((nxe+1)*(iq-1)/2)+ip;
                    num(2)=num(1)+1;
                    num(3)=floor((nxe+1)*(iq+1)/2)+ip;
                else
                    num(1)=(ip-1)*(nye+2)/2+(iq+1)/2;
                    num(2)=num(1)+(nye+2)/2;
                    num(3)=num(1)+1;
                end
                coord(1,1)=x_coords(ip);
                coord(1,2)=y_coords(floor((iq+1)/2));
                coord(2,1)=x_coords(ip+1)  ;
                coord(2,2)=y_coords((iq+1)/2);
                coord(3,1)=x_coords(ip)  ;
                coord(3,2)=y_coords((iq+3)/2);
            else
                if strcmp(dir,'x') || strcmp(dir,'r')
                    num(1)=floor((nxe+1)*iq/2)+ip+1;
                    num(2)=num(1)-1;
                    num(3)=floor((nxe+1)*(iq-2)/2)+ip+1;
                else
                    num(1)=floor(ip*(nye+2)/2)+floor((iq+2)/2);
                    num(2)=floor((ip-1)*(nye+2)/2)+floor((iq+1)/2)+1;
                    num(3)=num(1)-1;
                end
                coord(1,1)=x_coords(ip+1);
                coord(1,2)=y_coords(floor((iq+2)/2));
                coord(2,1)=x_coords(ip)   ;
                coord(2,2)=y_coords(floor((iq+2)/2));
                coord(3,1)=x_coords(ip+1) ;
                coord(3,2)=y_coords(floor(iq/2));
            end
        case 15
            if mod(iq,2)~=0
                if strcmp(dir,'x') || strcmp(dir,'r')
                    fac1=4*(4*nxe+1)*(iq-1)/2;
                    num(1)=fac1+4*ip-3;
                    num(2)=num(1)+1;
                    num(3)=num(1)+2;
                    num(4)=num(1)+3;
                    num(5)=num(1)+4;
                    num(6)=fac1+ 4*nxe+1+4*ip;
                    num(7)=fac1+ 8*nxe+1+4*ip;
                    num(8)=fac1+12*nxe+1+4*ip;
                    num(9)=fac1+16*nxe+1+4*ip;
                    num(10)=num(8)-1;
                    num(11)=num(7)-2;
                    num(12)=num(6)-3;
                    num(13)=num(12)+1;
                    num(14)=num(12)+2;
                    num(15)=num(11)+1;
                else
                    fac1=4*(2*nye+1)*(ip-1)+2*iq-1 ;
                    num(1)=fac1;
                    num(2)=fac1+2*nye+1;
                    num(3)=fac1+4*nye+2 ;
                    num(4)=fac1+6*nye+3 ;
                    num(5)=fac1+8*nye+4;
                    num(6)=fac1+6*nye+4 ;
                    num(7)=fac1+4*nye+4 ;
                    num(8)=fac1+2*nye+4;
                    num(9)=fac1+4 ;
                    num(10)=fac1+3 ;
                    num(11)=fac1+2 ;
                    num(12)=fac1+1;
                    num(13)=fac1+2*nye+2;
                    num(14)=fac1+4*nye+3;
                    num(15)=fac1+2*nye+3  ;
                end
                coord(1,1)=x_coords(ip);
                coord(1,2)=y_coords(floor((iq+1)/2));
                coord(5,1)=x_coords(ip+1)   ;
                coord(5,2)=y_coords(floor((iq+1)/2));
                coord(9,1)=x_coords(ip);
                coord(9,2)=y_coords(floor((iq+3)/2));
            else
                if strcmp(dir,'x') || strcmp(dir,'r')
                    fac1=4*(4*nxe+1)*(iq-2)/2;
                    num(1)=fac1+16*nxe+5+4*ip;
                    num(2)=num(1)-1;
                    num(3)=num(1)-2;
                    num(4)=num(1)-3;
                    num(5)=num(1)-4;
                    num(6)=fac1+12*nxe+1+4*ip;
                    num(7)=fac1+8*nxe+1+4*ip;
                    num(8)=fac1+4*nxe+1+4*ip;
                    num(9)=fac1+4*ip+1;
                    num(10)=num(8)+1;
                    num(11)=num(7)+2;
                    num(12)=num(6)+3;
                    num(13)=num(12)-1;
                    num(14)=num(12)-2;
                    num(15)=num(11)-1;
                else
                    fac1=4*(2*nye+1)*(ip-1)+2*iq+8*nye+5 ;
                    num(1)=fac1 ;
                    num(2)=fac1-2*nye-1;
                    num(3)=fac1-4*nye-2 ;
                    num(4)=fac1-6*nye-3 ;
                    num(5)=fac1-8*nye-4;
                    num(6)=fac1-6*nye-4 ;
                    num(7)=fac1-4*nye-4;
                    num(8)=fac1-2*nye-4;
                    num(9)=fac1-4;
                    num(10)=fac1-3 ;
                    num(11)=fac1-2 ;
                    num(12)=fac1-1;
                    num(13)=fac1-2*nye-2 ;
                    num(14)=fac1-4*nye-3;
                    num(15)=fac1-2*nye-3 ;
                end
                coord(1,1)=x_coords(ip+1);
                coord(1,2)=y_coords(floor((iq+2)/2));
                coord(5,1)=x_coords(ip)   ;
                coord(5,2)=y_coords(floor((iq+2)/2));
                coord(9,1)=x_coords(ip+1) ;
                coord(9,2)=y_coords(floor(iq/2));
            end
            coord(3,:)=pt5*(coord(1,:)+coord(5,:));
            coord(7,:)=pt5*(coord(5,:)+coord(9,:));
            coord(11,:)=pt5*(coord(9,:)+coord(1,:));
            coord(2,:)=pt5*(coord(1,:)+coord(3,:));
            coord(4,:)=pt5*(coord(3,:)+coord(5,:));
            coord(6,:)=pt5*(coord(5,:)+coord(7,:));
            coord(8,:)=pt5*(coord(7,:)+coord(9,:));
            coord(10,:)=pt5*(coord(9,:)+coord(11,:));
            coord(12,:)=pt5*(coord(11,:)+coord(1,:));
            coord(15,:)=pt5*(coord(7,:)+coord(11,:));
            coord(14,:)=pt5*(coord(3,:)+coord(7,:));
            coord(13,:)=pt5*(coord(2,:)+coord(15,:)) ;
        otherwise
            disp('Wrong number of nodes for triangular element');
    end
else
    nye=size(y_coords,1)-1;
    if strcmp(dir,'x') || strcmp(dir,'r')
        iq=floor((iel-1)/nxe)+1;
        ip=iel-(iq-1)*nxe;
    else
        ip=floor((iel-1)/nye)+1;
        iq=iel-(ip-1)*nye;
    end
    switch nod
        case 4
            if strcmp(dir,'x') || strcmp(dir,'r')
                num(1)=iq*(nxe+1)+ip;
                num(2)=(iq-1)*(nxe+1)+ip;
                num(3)=num(2)+1;
                num(4)=num(1)+1;
            else
                num(1)=(ip-1)*(nye+1)+iq+1;
                num(2)=num(1)-1;
                num(3)=ip*(nye+1)+iq;
                num(4)=num(3)+1;
            end
            coord(1:2,1)=x_coords(ip);
            coord(3:4,1)=x_coords(ip+1);
            coord(1,2)=y_coords(iq+1);
            coord(2:3,2)=y_coords(iq);
            coord(4,2)=coord(1,2);
        case 5
            if strcmp(dir,'x') || strcmp(dir,'r')
                num(1)=iq*(2*nxe+1)+ip;
                num(2)=(iq-1)*(2*nxe+1)+ip;
                num(3)=num(2)+1;
                num(4)=num(1)+1;
                num(5)=iq*(2*nxe+1)+ip-nxe;
            else
                num(1)=(ip-1)*(2*nye+1)+iq+1;
                num(2)=num(1)-1;
                num(3)=ip*(2*nye+1)+iq;
                num(4)=num(3)+1;
                num(5)=ip*(2*nye+1)+iq-nye;
            end
            coord(1:2,1)=x_coords(ip);
            coord(3:4,1)=x_coords(ip+1);
            coord(1,2)=y_coords(iq+1);
            coord(2:3,2)=y_coords(iq);
            coord(4,2)=coord(1,2);
            coord(5,:)=0.25*(coord(1,:)+coord(2,:)+coord(3,:)+coord(4,:));
        case 8
            if strcmp(dir,'x') || strcmp(dir,'r')
                num(1)=iq*(3*nxe+2)+2*ip-1 ;
                num(2)=iq*(3*nxe+2)+ip-nxe-1	;
                num(3)=(iq-1)*(3*nxe+2)+2*ip-1;
                num(4)=num(3)+1;
                num(5)=num(4)+1;
                num(6)=num(2)+1;
                num(7)=num(1)+2;
                num(8)=num(1)+1;
            else
                num(1)=(ip-1)*(3*nye+2)+2*iq+1;
                num(2)=num(1)-1;
                num(3)=num(1)-2;
                num(4)=(ip-1)*(3*nye+2)+2*nye+iq+1;
                num(5)=ip*(3*nye+2)+2*iq-1;
                num(6)=num(5)+1;
                num(7)=num(5)+2;
                num(8)=num(4)+1;
            end
            coord(1:3,1)=x_coords(ip);
            coord(5:7,1)=x_coords(ip+1);
            coord(4,1)=pt5*(coord(3,1)+coord(5,1));
            coord(8,1)=pt5*(coord(7,1)+coord(1,1));
            coord(1,2)=y_coords(iq+1);
            coord(7:8,2)=y_coords(iq+1);
            coord(3:5,2)=y_coords(iq);
            coord(2,2)=pt5*(coord(1,2)+coord(3,2));
            coord(6,2)=pt5*(coord(5,2)+coord(7,2));
        otherwise
            disp('Wrong number of nodes for quadrilateral element');
    end
end

