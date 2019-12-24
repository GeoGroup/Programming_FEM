function nf = emb_2d_bc(nx1,nx2,ny1,ny2)
% generates the nf array for a 2-d slope geometry considering boundary condition.
nye=ny1+ny2;
nn=(3*nye+2)*nx1+2*nye+1+(3*ny2+2)*nx2;
nodof=2;
nf=zeros(nodof,nn);
nm=0;
ic=0;
for i=1:2*nye
    nm=nm+1;
    nf(1,nm)=0;
    ic=ic+1;
    nf(2,nm)=ic;
end
nm=nm+1;
nf(1,nm)=0;
nf(2,nm)=0;
for j=1:nx1-1
    for i=1:nye
        nm=nm+1;
        ic=ic+1;
        nf(1,nm)=ic;
        ic=ic+1;
        nf(2,nm)=ic;
    end
    nm=nm+1;
    nf(1,nm)=0;
    nf(2,nm)=0;
    for i=1:2*nye
        nm=nm+1;
        ic=ic+1;
        nf(1,nm)=ic;
        ic=ic+1;
        nf(2,nm)=ic;
    end
    nm=nm+1;
    nf(1,nm)=0;
    nf(2,nm)=0;
end
for i=1:nye
    nm=nm+1;
    ic=ic+1;
    nf(1,nm)=ic;
    ic=ic+1;
    nf(2,nm)=ic;
end
nm=nm+1;
nf(1,nm)=0;
nf(2,nm)=0;
for i=1:2*ny1
    nm=nm+1;
    ic=ic+1;
    nf(1,nm)=ic;
    ic=ic+1;
    nf(2,nm)=ic;
end
if (nx2==0)
    for i=1:2*ny2
        nm=nm+1;
        nf(1,nm)=0;
        ic=ic+1;
        nf(2,nm)=ic;
    end
    nm=nm+1;
    nf(1,nm)=0;
    nf(2,nm)=0;
else
    for i=1:2*ny2
        nm=nm+1;
        ic=ic+1;
        nf(1,nm)=ic;
        ic=ic+1;
        nf(2,nm)=ic;
    end
    nm=nm+1;
    nf(1,nm)=0;
    nf(2,nm)=0;
    for j=1:nx2
        for i=1:ny2
            nm=nm+1;
            ic=ic+1;
            nf(1,nm)=ic;
            ic=ic+1;
            nf(2,nm)=ic;
        end
        nm=nm+1;
        nf(1,nm)=0;
        nf(2,nm)=0;
        for i=1:2*ny2
            nm=nm+1;
            if(j==nx2)
                nf(1,nm)=0;
            else
                ic=ic+1;
                nf(1,nm)=ic;
            end
            ic=ic+1;
            nf(2,nm)=ic;
        end
        nm=nm+1;
        nf(1,nm)=0;
        nf(2,nm)=0;
    end
end
end

