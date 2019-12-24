function f = mocouf(phi,c,sigm,dsbar,theta)
% ! This subroutine calculates the value of the yield function
% ! for a Mohr-Coulomb material (phi in degrees).
one=1.0;d3=3.0;d4=4.0;
d180=180.0;
phir=phi*d4*atan(one)/d180;
snph=sin(phir) ;
csph=cos(phir) ;
csth=cos(theta);
snth=sin(theta);
f=snph*sigm+dsbar*(csth/sqrt(d3)-snth*snph/d3)-c*csph;

end

