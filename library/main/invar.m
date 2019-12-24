function [sigm,dsbar,theta]=invar(stress)
% ! This subroutine forms the stress invariants in 2- or 3-d.
zero=0.0;small=1.e-10;one=1.0;two=2.0;
three=3.0;six=6.0;thpt5=13.5;
nst=size(stress,1);
switch nst
    case 4
        sx=stress(1);
        sy=stress(2);
        txy=stress(3);
        sz=stress(4);
        sigm=(sx+sy+sz)/three;
        dsbar=sqrt((sx-sy)^2+(sy-sz)^2+(sz-sx)^2+six*txy^2)/sqrt(two);
        if dsbar<small
            theta=zero;
        else
            dx=(two*sx-sy-sz)/three;
            dy=(two*sy-sz-sx)/three;
            dz=(two*sz-sx-sy)/three;
            xj3=dx*dy*dz-dz*txy^2;
            sine=-thpt5*xj3/dsbar^3;
            if sine>=one;sine=one;end
            if sine<-one;sine=-one;end
            theta=asin(sine)/three ;
        end
    case 6
        sq3=sqrt(three);
        s1=stress(1);
        s2=stress(2);
        s3=stress(3);
        s4=stress(4);
        s5=stress(5);
        s6=stress(6);
        sigm=(s1+s2+s3)/three;
        d2=((s1-s2)^2+(s2-s3)^2+(s3-s1)^2)/six+s4*s4+s5*s5+s6*s6;
        ds1=s1-sigm ;
        ds2=s2-sigm  ;
        ds3=s3-sigm;
        d3=ds1*ds2*ds3-ds1*s5*s5-ds2*s6*s6-ds3*s4*s4+two*s4*s5*s6;
        dsbar=sq3*SQRT(d2);
        if dsbar<small
            theta=zero;
        else
            sine=-three*sq3*d3/(two*SQRT(d2)^3);
            if sine>=one;sine=one;end
            if sine<-one;sine=-one;end
            theta=asin(sine)/three ;
        end
    otherwise
        disp('wrong size for nst in invar');
end

end

