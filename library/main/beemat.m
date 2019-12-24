function bee =beemat( deriv,ih )
% This subroutine forms the bee matrix in 2-d (ih=3 or 4) or 3-d (ih=6).
nod=size(deriv,2);

switch ih
    case {3,4}
        bee=zeros(ih,nod*2);
        for m=1:nod
            k=2*m;
            l=k-1;
            x=deriv(1,m);
            y=deriv(2,m);
            bee(1,l)=x;
            bee(3,k)=x;
            bee(2,k)=y;
            bee(3,l)=y;
        end
    case 6
        bee=zeros(ih,nod*3);
           for m=1:nod
               n=3*m;
               k=n-1;
               l=k-1;
               x=deriv(1,m);
               y=deriv(2,m);
               z=deriv(3,m);
               bee(1,l)=x;
               bee(4,k)=x;
               bee(6,n)=x;
               bee(2,k)=y;
               bee(4,l)=y;
               bee(5,n)=y;
               bee(3,n)=z;
               bee(5,k)=z;
               bee(6,l)=z;
           end
    otherwise
        disp('wrong dimension for nst in bee matrix');
end

end

