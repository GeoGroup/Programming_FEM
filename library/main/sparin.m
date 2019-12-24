function kv = sparin( kv,kdiag )
% Cholesky factorisation on a symmetric  skyline global matrix
 n=size(kdiag,1)  ;
 kv(1)=sqrt(kv(1));
 
  for i=2:n
   ki=kdiag(i)-i;
   l=kdiag(i-1)-ki+1;
   for j=l:i
     x=kv(ki+j);
     kj=kdiag(j)-j;
     if j~=1
       ll=kdiag(j-1)-kj+1;
       ll=max(l,ll);
       if (ll~=j)
         m=j-1;
         for k=ll:m 
           x=x-kv(ki+k)*kv(kj+k) ;
         end
       end
     end
     kv(ki+j)=x/kv(kj+j);
   end
   kv(ki+i)=sqrt(x);
  end
end

