function [bee,radius] = bmat_nonaxi(coord,deriv,fun,iflag,lth)
nod=size(deriv,2);
bee=zeros(6,3*nod);
radius=fun'*coord(:,1);
for m=1:nod
    n=3*m;
    k=n-1;
    l=k-1;
    bee(1,l)=deriv(1,m);
    bee(2,k)=deriv(2,m);
    bee(3,l)=fun(m)/radius;
    bee(3,n)=iflag*lth*bee(3,l) ;
    bee(4,l)=deriv(2,m);
    bee(4,k)=deriv(1,m);
    bee(5,k)=-iflag*lth*fun(m)/radius ;
    bee(5,n)=deriv(2,m);
    bee(6,l)=bee(5,k) ;
    bee(6,n)=deriv(1,m)-fun(m)/radius;
end
end

