function [dq1,dq2,dq3]=mocouq(psi,dsbar,theta)
% ! This subroutine forms the derivatives of a Mohr-Coulomb potential
% ! function with respect to the three invariants (psi in degrees).
pt49=0.49;pt5=0.5;one=1.0;d3=3.0;d4=4.0;zero=0;
d180=180.0;
psir=psi*d4*atan(one)/d180 ;
snth=sin(theta) ;
snps=sin(psir);
sq3=sqrt(d3)  ;
dq1=snps;
if abs(snth) > pt49
    c1=one;
    if snth  < zero
        c1=-one;
    end
    dq2=(sq3*pt5-c1*snps*pt5/sq3)*sq3*pt5/dsbar ;
    dq3=zero;
else
    csth=cos(theta);
    cs3th=cos(d3*theta);
    tn3th=tan(d3*theta);
    tnth=snth/csth;
    dq2=sq3*csth/dsbar*((one+tnth*tn3th)+snps*(tn3th-tnth)/sq3)*pt5;
    dq3=pt5*d3*(sq3*snth+snps*csth)/(cs3th*dsbar*dsbar);
end
end

