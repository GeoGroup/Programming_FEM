function nf = bc_rect(nxe,nye,dir)
% ! This subroutine generates the nf array for a rectangular mesh
% ! of 8-node quadrilaterals fully fixed on the base and with
% ! vertical rollers on the left and right sides. Nodes numbered
% ! in the x- or y-direction
nf=0;
if strcmp(dir,'y')
    nm=0;
    ic=0;
    for j=1:2*nye
        nm=nm+1;
        nf(1,nm)=0;
        ic=ic+1;
        nf(2,nm)=ic;
    end
    nm=nm+1;
    nf(1,nm)=0;
    nf(2,nm)=0;
    %   !
    for i=1:nxe-1
        for j=1:nye
            nm=nm+1;
            ic=ic+1;
            nf(1,nm)=ic;
            ic=ic+1;
            nf(2,nm)=ic;
        end
        nm=nm+1;
        nf(1,nm)=0;
        nf(2,nm)=0 ;
        %   !
        for j=1:2*nye;
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
    %   !
    for j=1:nye
        nm=nm+1;
        ic=ic+1;
        nf(1,nm)=ic;
        ic=ic+1;
        nf(2,nm)=ic;
    end
    nm=nm+1;
    nf(1,nm)=0;
    nf(2,nm)=0 ;
    %   !
    for j=1:2*nye
        nm=nm+1;
        nf(1,nm)=0;
        ic=ic+1;
        nf(2,nm)=ic;
    end
    nm=nm+1;
    nf(1,nm)=0;
    nf(2,nm)=0;
else
    nm=0;
    ic=0;
    for j=1:nye
        nm=nm+1;
        nf(1,nm)=0;
        ic=ic+1;
        nf(2,nm)=ic;
        for i=1:2*nxe-1
            nm=nm+1;
            ic=ic+1;
            nf(1,nm)=ic;
            ic=ic+1;
            nf(2,nm)=ic;
        end
        nm=nm+1;
        nf(1,nm)=0;
        ic=ic+1;
        nf(2,nm)=ic;
        % !
        nm=nm+1;
        nf(1,nm)=0;
        ic=ic+1;
        nf(2,nm)=ic;
        for i=1:nxe-1;
            nm=nm+1;
            ic=ic+1;
            nf(1,nm)=ic;
            ic=ic+1;
            nf(2,nm)=ic;
        end
        nm=nm+1;
        nf(1,nm)=0;
        ic=ic+1;
        nf(2,nm)=ic;
    end
    % !
    for i=1:2*nxe+1
        nm=nm+1;
        nf(1,nm)=0;
        nf(2,nm)=0;
    end
end
end

