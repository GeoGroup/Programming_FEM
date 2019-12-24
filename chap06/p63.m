% 2D Viscoplastic 
clear;clc;
addpath ../geom
addpath ../main
zero=0;two=2; dt=1.0e15;d3=3.0;d4=4.0;d6=6;one=1.0;%nprops=6;
penalty=1e20;start_dt=1.e15;
ndim=2;nodof=2;nip=4;nod=8;nst=4;
ndof=nodof*nod;
global nxe nye;
element='quadrilateral';
[filename, pathname] = uigetfile( ...
       {'*.dat';'*.txt';'*.*'}, ...
        'Pick a file');
name=[pathname,filename];
%     [type_2d,element,nod,dir,nxe,nye,nip
fid=fopen(name,'rt');
tmp=sscanf(fgetl(fid),'%d %d %d');
nxe=tmp(1);
nye=tmp(2);
nbo2=tmp(3);
[ nels,nn ] = mesh_size( element,nod );
np_types=sscanf(fgetl(fid),'%d');
prop=sscanf(fgetl(fid),'%f %f %f %f %f %f');


etype=ones(nels,1);
if np_types>1
    for i=1:nels
        etype(i)=sscanf(fgetl(fid),'%d');
    end
end
qs=sscanf(fgetl(fid),'%f');
nf=bc_rect(nxe,nye,'y');
neq=max(nf(:));

str='x_coords=sscanf(fgetl(fid),''';
for i=1:nxe
    str=[str,'%f '];
end
str=[str,''');'];
eval(str);
str='y_coords=sscanf(fgetl(fid),''';
for i=1:nye
    str=[str,'%f '];
end
str=[str,''');'];
eval(str);

kdiag=zeros(neq,1);
g_num=zeros(nod,nels);
g_coord=zeros(ndim,nn);
g_g=zeros(ndof,nels);
g=zeros(ndof,1); %每个单元的自由度序号
% !-----------------------loop the elements to find global arrays sizes-----
for iel=1:nels %element_1
    [ num,coord ] = geom_rect(iel,element,x_coords,y_coords,nod ,'y');
    [ g ] = num_to_g(num,nf);
    g_num(:,iel)=num;
    g_coord(:,num)=coord';
    g_g(:,iel)=g ;
    kdiag = fkdiag(kdiag, g);
end
for i=2:neq
    kdiag(i)=kdiag(i)+kdiag(i-1);
end
disp(['There are ',num2str(neq) ' equations and the skyline storage is ',num2str(kdiag(neq))]);
dt=start_dt;
for i=1:np_types
   phi=prop(1,i);
   snph=sin(phi*pi/180);
   e=prop(5,i);
   v=prop(6,i) ;
   ddt=d4*(one+v)*(one-two*v)/(e*(one-two*v+snph^2)); %计算dt
   if(ddt<dt);dt=ddt;end
end 
[points, weights]= sample(element,nip );
kv=zeros(kdiag(neq),1);    
gravlo=zeros(neq,1);
% % !-----------------------element stiffness integration and assembly--------
for iel=1:nels % elements_2
    dee=deemat(prop(5,etype(iel)),prop(6,etype(iel)),nst);
    num=g_num(:,iel);
    coord=g_coord(:,num)';
    g=g_g(:,iel);
    km=zeros(ndof,ndof);
    eld=zeros(ndof,1);
    % gauss integration
    for i=1:nip
        fun=shape_fun(points,i,ndim,nod);
        der=shape_der(points,i,ndim,nod);
        jac=der*coord;
        Det=det(jac);
        deriv=jac\der;
        bee=beemat( deriv,nst );
        km=km+bee'*dee*bee*Det*weights(i);
        eld(nodof:nodof:ndof)=eld(nodof:nodof:ndof)+fun(:)*Det*weights(i);
    end
    kv=fsparv(kv,km,g,kdiag);
    for i=1:ndof
        if g(i)~=0
            gravlo(g(i))=gravlo(g(i))-eld(i)*prop(4,etype(iel));
        end
    end
end % elements_2
kvc=kv;
% !-----------------------surcharge loads-----------------------------------
for i=1:nxe
   i3=g_num(3,(i-1)*nye+1);
   i4=g_num(4,(i-1)*nye+1);
   i5=g_num(5,(i-1)*nye+1);
   qq=(x_coords(i+1)-x_coords(i))*qs;
   gravlo(nf(2,i3))=gravlo(nf(2,i3))-qq/d6;
   gravlo(nf(2,i4))=gravlo(nf(2,i4))-qq*two/d3;
   gravlo(nf(2,i5))=gravlo(nf(2,i5))-qq/d6;
end 
% !-----------------------factorise equations-------------------------------
kv=sparin(kv,kdiag);
gravlo = spabac(kv,gravlo,kdiag);
% !-----------------------set up initial stresses---------------------------
tensor=zeros(nst,nip,nels);
for iel=1:nels %elements_3
    dee=deemat(prop(5,etype(iel)),prop(6,etype(iel)),nst);
    g=g_g(:,iel);
    for i=1:length(g)
        if g(i)~=0
            eld(i)=gravlo(g(i));
        else
            eld(i)=0;
        end
    end
    num=g_num(:,iel);
    coord=(g_coord(:,num))';
    for i=1:nip %int_pts_2
        der=shape_der(points,i,ndim,nod);
        jac=der*coord;
        Det=det(jac);
        deriv=jac\der;
        bee=beemat( deriv,nst );
        sigma=dee*(bee*eld);
        tensor(:,i,iel)=sigma;
    end       %int_pts_2
end %elements_3
% !-----------------------fixed displacement data and factorise equations----
fixed_freedoms=2*nbo2+1;
node=zeros(fixed_freedoms,1);
no=zeros(fixed_freedoms,1);
% storkv=zeros(fixed_freedoms,1);
node(1)=1;k=1;
for i=1:nbo2
    k=k+2*nye+1;
    node(2*i)=k;
    k=k+nye+1;
    node(2*i+1)=k;
end
for i=1:fixed_freedoms
   no(i)=nf(2,node(i));
end
kv=kvc;
kv(kdiag(no))=kv(kdiag(no))+penalty;
storkv=kv(kdiag(no));
kv=sparin(kv,kdiag);
% !-----------------------load increment loop-------------------------------
tmp=sscanf(fgetl(fid),'%f %d');
tol=tmp(1);limit=tmp(2);
tmp=sscanf(fgetl(fid),'%d %f');
incs=tmp(1);presc=tmp(2);
oldis=zeros(neq,1);totd=zeros(neq,1);
devp=zeros(nst,1);
for iy=1:incs %incs %disp_incs
    bdylds=zeros(neq,1);react=zeros(neq,1);
    iters=0;
    evpt=zeros(nst,nip,nels);
    % !-----------------------plastic iteration loop----------------------------
    while 1 %its
% for iii=1:1
        iters=iters+1;loads=zeros(neq,1);
        loads=loads+bdylds;
        for i=1:fixed_freedoms
            loads(no(i))=storkv(i)*presc;
        end
        loads = spabac(kv,loads,kdiag);
        % !-----------------------check plastic convergence-------------------------
        [converged,oldis,res] = checon(loads,oldis,tol);
        
%         disp(['res equal to ',num2str(res)]);

        if iters==1
            converged=0;
        end
        if converged==1 || iters==limit
            bdylds=zeros(neq,1);
        end
        % !-----------------------go round the Gauss Points ------------------------
        for iel=1:nels %elements_4
            phi=prop(1,etype(iel));
            c=prop(2,etype(iel));
            psi=prop(3,etype(iel)) ;
            dee=deemat(prop(5,etype(iel)),prop(6,etype(iel)),nst);
            bload=zeros(ndof,1);rload=zeros(ndof,1);
            num=g_num(:,iel);
            coord=g_coord(:,num)';
            g=g_g(:,iel);
            eld=zeros(length(g),1);
            for ii=1:length(g) %eld=loads(g) 由于matlab下表从1 0自由度比较麻烦
                if g(ii)~=0
                    eld(ii)=loads(g(ii));
                else
                    eld(ii)=0;
                end
            end
            for i=1:nip %gauss_pts_4
                der=shape_der(points,i,ndim,nod);
                jac=der*coord;
                Det=det(jac);
                deriv=jac\der;
                bee=beemat( deriv,nst );
                eps=bee*eld;
                eps=eps-evpt(:,i,iel);
                sigma=dee*eps;
                stress=sigma+tensor(:,i,iel) ;
                [sigm,dsbar,lode_theta]=invar(stress);
                %  !-----------------------check whether yield is violated-------------------
                f=mocouf(phi,c,sigm,dsbar,lode_theta);
                if converged==1 && iters==limit
                    devp=stress;
                else
                    if f>=zero
                        [dq1,dq2,dq3]=mocouq(psi,dsbar,lode_theta);
                        [m1,m2,m3]= formm(stress);
                        flow=f*(m1*dq1+m2*dq2+m3*dq3);
                        erate=flow*stress;
                        evp=erate*dt;
                        evpt(:,i,iel)=evpt(:,i,iel)+evp;
                        devp=dee*evp;
                    end
                end
                if f>=0 || (converged==1 || iters==limit)
                    eload=bee'*devp;
                    bload=bload+eload*Det*weights(i);
                end
                % !-----------------------update the Gauss Point stresses-------------------
                if converged ==1 || iters==limit
                    tensor(:,i,iel)=stress;
                    rload=rload+bee'*stress*Det*weights(i);
                end
            end %gauss_pts_4
            % !-----------------------compute the total bodyloads vector----------------
            for jj=1:length(g)
                if g(jj)~=0
                    bdylds(g(jj))=bdylds(g(jj))+bload(jj);
                    react(g(jj))=react(g(jj))+rload(jj);
                end
            end
        end %elements_4
        if converged==1 || iters==limit
            break;
        end
% %         disp([num2str(iters),'   iters   ',num2str(res)]);
    end %its
    
    totd=totd+loads;
    pr=zero;
    for i=1:fixed_freedoms
        pr=pr+react(no(i));
    end
    pr=pr/x_coords(nbo2+1);
    pav=zero;
    for i=1:nbo2
        pav=pav+tensor(2,1,(i-1)*nye+1)+tensor(2,2,(i-1)*nye+1);
    end
    pav=pav/(two*nbo2);
    disp([num2str(iy),'   ',num2str(-totd(1)),'  ',num2str(-pr),'  ',num2str(-pav),'   ',num2str(iters)]);
    if iters==limit
        break;
    end
    
%     fprintf('%2d %3d %8.5g %3d \n',iy,qinc(iy),loads(nf(2,node(1))),iters);
% %     disp([num2str(iy),num2str(ptot),num2str(totd(nf(2,node(1)))),num2str(iters)]);
end %disp_incs
