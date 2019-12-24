function kv = fsparv(kv,km,g,kdiag)
% Assembles element matrices into a symmetric skyline
idof=length(g);
for i=1:idof
    k=g(i);
    if k~=0
        for j=1:idof
            if g(j)~=0
                iw=k-g(j);
                if iw>=0
                    ival=kdiag(k)-iw;
                    kv(ival)=kv(ival)+km(i,j) ;
                end
            end
        end
    end
end
end

