function loads = spabac(kv,loads,kdiag)
% This subroutine performs Cholesky forward and back-substitution
% on a symmetric skyline global matrix.
 n=size(kdiag,1);
 loads(1)=loads(1)/kv(1);
 for i=2:n
   ki=kdiag(i)-i;
   l=kdiag(i-1)-ki+1 ;
   x=loads(i);
   if (l~=i)
     m=i-1;
     for j=l:m 
       x=x-kv(ki+j)*loads(j);
     end
   end
   loads(i)=x/kv(ki+i);
 end
 for it=2:n
   i=n+2-it;
   ki=kdiag(i)-i;
   x=loads(i)/kv(ki+i);
   loads(i)=x;
   l=kdiag(i-1)-ki+1;
   if (l~=i)
     m=i-1;
     for k=l:m
       loads(k)=loads(k)-x*kv(ki+k);
     end
   end
 end
 loads(1)=loads(1)/kv(1);


end

