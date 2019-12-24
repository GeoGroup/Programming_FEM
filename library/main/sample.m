function [points, weights]= sample(element,nip )
%  This subroutine returns the local coordinates and weighting coefficients
% of the integrating points.
points=0;weights=8;
root3=1/sqrt(3);r15=0.2*sqrt(15.0);
w=[5.0/9.0;8.0/9.0;5.0/9.0];
v=[5.0/9.0*w;8.0/9.0*w;5.0/9.0*w];
switch element
    case 'triangle'
        switch nip
            case 1
                points=[1/3,1/3];weights=0.5;
            case 3
                points=[0.5 0.5
                    0.5 0
                    0 0.5];
                weights=ones(nip,1)/6;
            case 4
                points=[0.6 0.2
                    0.2 0.6
                    0.2 0.2
                    1/3 1/3];
                weights=[0.520833333333333;0.520833333333333;0.520833333333333;-0.5625]/2;
            case 6
                points=[0.816847572980459 0.091576213509771
                    0.091576213509771 0.816847572980459
                    0.091576213509771 0.091576213509771
                    0.108103018168070 0.445948490915965
                    0.445948490915965 0.108103018168070
                    0.445948490915965 0.445948490915965];
                weights=[0.109951743655322
                    0.109951743655322
                    0.109951743655322
                    0.223381589678011
                    0.223381589678011
                    0.223381589678011]/2;
            case 12
                s=zeros(12,2);
                s(1,1)= 0.873821971016996;
                s(1,2)= 0.063089014491502;
                s(2,1)= 0.063089014491502;
                s(2,2)= 0.873821971016996;
                s(3,1)= 0.063089014491502;
                s(3,2)= 0.063089014491502;
                s(4,1)= 0.501426509658179;
                s(4,2)= 0.249286745170910;
                s(5,1)= 0.249286745170910;
                s(5,2)= 0.501426509658179;
                s(6,1)= 0.249286745170910;
                s(6,2)= 0.249286745170910;
                s(7,1) =0.053145049844817;
                s(7,2) =0.310352451033784;
                s(8,1) =0.310352451033784;
                s(8,2) =0.053145049844817;
                s(9,1) =0.053145049844817;
                s(9,2) =0.636502499121398;
                s(10,1)=0.310352451033784;
                s(10,2)=0.636502499121398;
                s(11,1)=0.636502499121398;
                s(11,2)=0.053145049844817;
                s(12,1)=0.636502499121398;
                s(12,2)=0.310352451033784;
                points=s;
                wt=zeros(12,1);
                wt(1:3)=0.050844906370207;
                wt(4:6)=0.116786275726379;
                wt(7:12)=0.082851075618374;
                wt=0.5*wt;
                weights=wt;
            otherwise
                disp('wrong number of integrating points for a triangle');
        end
    case 'quadrilateral'
        switch nip
            case 1
                points=[0 0];
                weights=4;
            case 4
                points=[-root3 root3
                    root3 root3
                    -root3 -root3
                    root3 -root3];
                weights=[1;1;1;1];
            case 9
                s=zeros(9,2);
                s(1:3:7,1)=-r15;
                s(2:3:8,1)=0.0;
                s(3:3:9,1)=r15;
                s(1:3,2)  =r15;
                s(4:6,2)  =0.0;
                s(7:9,2)  =-r15;
                wt= v;
                points=s;
                weights=wt;
            otherwise
                disp('wrong number of integrating points for a quadrilateral');
        end  
    case 'tetrahedron'
    switch nip
        case 1
            points=[0.25 0.25 0.25];
            weights=1/6;
        case 4
            s=zeros(4,3);
            s(1,1)=0.58541020;
            s(1,2)=0.13819660;
            s(1,3)=s(1,2);
            s(2,2)=s(1,1);
            s(2,3)=s(1,2);
            s(2,1)=s(1,2);
            s(3,3)=s(1,1);
            s(3,1)=s(1,2);
            s(3,2)=s(1,2);
            s(4,1)=s(1,2);
            s(4,2)=s(1,2);
            s(4,3)=s(1,2);
            wt=0.25/6.0*ones(4,1);
            points=s;
            weights=wt;
        otherwise
            disp('wrong number of integrating points for a tetrahedron');
    end
    case 'hexahedron'  
    switch nip
        case 1
            points=[0 0 0];
            weights=8;
        case 8
            points=[root3 root3 root3
                root3 root3 -root3
                root3 -root3 root3
                root3 -root3 -root3
                -root3 root3 root3
                -root3 -root3 root3
                -root3 root3 -root3
                -root3 -root3 -root3];
            weights=ones(nip,1);  
        case 14
            b=0.795822426;
            c=0.758786911;
            wt=zeros(14,1);
            wt(1:6)=0.886426593;
            wt(7:14)=0.335180055;
            weights=wt;
            s=zeros(14,3);
            s(1,1)=-b;
            s(2,1)=b;
            s(3,2)=-b;
            s(4,2)=b;
            s(5,3)=-b;
            s(6,3)=b;
            s(7:14,:)=c;
            s(7,1)=-c;
            s(7,2)=-c;
            s(7,3)=-c;
            s(8,2)=-c;
            s(8,3)=-c;
            s(9,1)=-c;
            s(9,3)=-c;
            s(10,3)=-c;
            s(11,1)=-c;
            s(11,2)=-c;
            s(12,2)=-c;
            s(13,1)=-c;
            points=s;
        case 27
            weights=[5.0/9.0*v;8.0/9.0*v;5.0/9.0*v];
            s(1:3:7,1)=-r15;
            s(2:3:8,1)=0.0;
            s(3:3:9,1)=r15;
            s(1:3,3)=r15;
            s(4:6,3)=0.0;
            s(7:9,3)=-r15;
            s(1:9,2)=-r15;
            s(10:3:16,1)=-r15;
            s(11:3:17,1)=0.0;
            s(12:3:18,1)=r15;
            s(10:12,3)=r15;
            s(13:15,3)=0.0;
            s(16:18,3)=-r15;
            s(10:18,2)=0.0;
            s(19:3:25,1)=-r15;
            s(20:3:26,1)=0.0;
            s(21:3:27,1)=r15;
            s(19:21,3)= r15;
            s(22:24,3)=0.0;
            s(25:27,3)=-r15;
            s(19:27,2)= r15;
            points=s;
        otherwise
            disp('wrong number of integrating points for a hexahedron');
    end
end

end

