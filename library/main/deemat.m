function dee = deemat(e,v,ih)
%! This subroutine returns the elastic dee matrix for ih=3 (plane strain),
%! ih=4 (axisymmetry or plane strain elastoplasticity) or ih=6
%! (three dimensions).
v1=1-v;
c=e/((1+v)*(1-2*v));
switch ih
    case 3
        dee=[v1*c v*c 0
            v*c v1*c 0
            0 0 0.5*c*(1-2*v)];
    case 4
        dee=[v1*c  v*c        0          v*c
            v*c    v1*c       0          v*c
            0       0    0.5*c*(1-2*v)    0
            v*c     v*c       0           v1*c];
        
    case 6
        v2=v/(1-v);
        vv=(1-2*v)/(1-v)*0.5;
        dee=[1  v2  v2  0  0   0
            v2  1   v2  0  0   0
            v2  v2  1   0  0   0
            0   0   0  vv  0   0
            0   0   0   0  vv  0
            0   0   0   0  0   vv]...
            *e/(2*(1+v)*vv);
end
end

