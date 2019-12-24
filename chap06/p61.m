% 2D Viscoplastic 
clear;clc;
addpath ../geom
addpath ../main
zero=0;two=2;dt=1.0e15;d3=3.0;d4=4.0;one=1.0;
penalty=1e20;
ndim=2;nodof=2;nip=4;nod=8;nst=4;
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
np_types=tmp(3);
[ nels,nn ] = mesh_size( element,nod );
etype=ones(nels,1);
if np_types>1
    for i=1:nels
        etype(i)=sscanf(fgetl(fid),'%d');
    end
end
prop=sscanf(fgetl(fid),'%f %f %f');

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
 
nf=ones(nodof,nn);
nr=sscanf(fgetl(fid),'%d');
for i=1:nr
    tmp=sscanf(fgetl(fid),'%d %d %d');
    nf(:,tmp(1))=tmp(2:3);
end
nf=formnf(nf);
neq=max(nf(:));

ndof=nod*nodof;
kdiag=zeros(neq,1);
g_num=zeros(nod,nels);
g_coord=zeros(ndim,nn);
g_g=zeros(ndof,nels);
g=zeros(ndof,1); %每个单元的自由度序号

for iel=1:nels
    [ num,coord ] = geom_rect(iel,element,x_coords,y_coords,nod ,dir);
    [ g ] = num_to_g(num,nf);
    g_num(:,iel)=num;
    g_coord(:,num)=coord';
    g_g(:,iel)=g ;
    kdiag = fkdiag(kdiag, g);
end
for i=2:neq
    kdiag(i)=kdiag(i)+kdiag(i-1);
end
disp(['There are ',num2str(neq) ' ,equations and the skyline storage is ',num2str(kdiag(neq))]);
% !-----------------------element stiffness integration and assembly--------
[points, weights]= sample(element,nip );
ih=nst;%ih=3 (plane strain)
%ih=4 (axisymmetry or plane strain elastoplasticity) 
% ih=6(three dimensions)
etype=ones(nels,1);
gc=1;
kv=zeros(kdiag(neq),1);
for iel=1:nels % elements_2
    ddt=d4*(one+prop(3,etype(iel)))/(d3*prop(2,etype(iel)));
    if ddt<dt
        dt=ddt;
    end
    dee=deemat(prop(2,etype(iel)),prop(3,etype(iel)),ih);
    num=g_num(:,iel);
    coord=g_coord(:,num)';
    g=g_g(:,iel);
    km=zeros(ndof,ndof);
    % gauss integration
    for i=1:nip
        fun=shape_fun(points,i,ndim,nod);
        der=shape_der(points,i,ndim,nod);
        jac=der*coord;
        Det=det(jac);
        deriv=jac\der;
%         jac=inv(jac);
%         deriv=jac*der;
        bee=beemat( deriv,ih );
        km=km+bee'*dee*bee*Det*weights(i)*gc(1);
    end
    kv=fsparv(kv,km,g,kdiag);
end % elements_2
% !-----------------------read load weightings and factorise equations------
loaded_nodes=sscanf(fgetl(fid),'%d');
if loaded_nodes>0
    node=zeros(loaded_nodes,1);
    value=zeros(loaded_nodes,ndim);
    for i=1:loaded_nodes
        tmp=sscanf(fgetl(fid),'%d %f %f');
        node(i)=tmp(1);
        value(i,:)=tmp(2:3);
    end
end
kv=sparin(kv,kdiag);
% !-----------------------load increment loop-------------------------------
loads=zeros(neq,1);bdylds=loads;oldis=loads;totd=loads;
tmp=sscanf(fgetl(fid),'%f %d');
tol=tmp(1);limit=tmp(2);
incs=sscanf(fgetl(fid),'%d');
qinc=zeros(incs,1);
for i=1:incs
    qinc(i)=sscanf(fgetl(fid),'%f');
end
tensor=zeros(nst,nip,nels);
ptot=0;devp=zeros(nst,1);
disp(['增量步 ','荷载增量  ','位移增量  ','迭代次数']);
for iy=1:incs
    ptot=ptot+qinc(iy);
    iters=0;bdylds=zeros(neq,1);evpt=zeros(nst,nip,nels);
    % !-----------------------plastic iteration loop----------------------------
    while 1 %its
% for kk=1:250
        iters=iters+1;loads=zeros(neq,1);
        for i=1:loaded_nodes
            index=nf(:,node(i));
            for j=1:ndim
                if index(j)~=0
                    loads(index(j))=value(i,j)*qinc(iy);
                end
            end
        end
        loads=loads+bdylds;
        loads = spabac(kv,loads,kdiag);
        % !-----------------------check plastic convergence-------------------------
        [converged,oldis,res] = checon(loads,oldis,tol);
        if iters==1
            converged=0;
        end
        if converged==1 || iters==limit
            bdylds=zeros(neq,1);
        end
        % !-----------------------go round the Gauss Points ------------------------
        for iel=1:nels %elements_3
            dee=deemat(prop(2,etype(iel)),prop(3,etype(iel)),ih);
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
            bload=zeros(ndof,1);
            for i=1:nip %gauss_pts_2
                der=shape_der(points,i,ndim,nod);
                jac=der*coord;
                Det=det(jac);
                deriv=jac\der;
                bee=beemat( deriv,ih );
                eps=bee*eld;
                eps=eps-evpt(:,i,iel);
                sigma=dee*eps;
                stress=sigma+tensor(:,i,iel) ;
                [sigm,dsbar,theta]=invar(stress);
                %  !-----------------------check whether yield is violated-------------------
                f=dsbar-sqrt(d3)*prop(1,etype(iel));
                if converged==1 && iters==limit
                    devp=stress;
                else
                    if f>=zero
                        dq1=zero;
                        dq2=d3/two/dsbar;
                        dq3=zero;
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
                end
            end %gauss_pts_2
% !-----------------------compute the total bodyloads vector----------------
            for jj=1:length(g)
                if g(jj)~=0
                    bdylds(g(jj))=bdylds(g(jj))+bload(jj);
                end
            end
        end %elements_3
        if converged==1 || iters==limit
            break;
        end
%         disp([num2str(iters),'   iters   ',num2str(res)]);
% end
    end %its
    totd=totd+loads;
    fprintf('%2d %8d %12.5g %6d \n',iy,qinc(iy),loads(nf(2,node(1))),iters);
%     disp([num2str(iy),num2str(ptot),num2str(totd(nf(2,node(1)))),num2str(iters)]);
end %load_incs

