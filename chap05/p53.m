clear;clc;
addpath ../geom
addpath ../main
penalty=1e20;
ndim=3;nodof=3;
global nxe nye nze;
nst=6; %The row number of B Matrix
element='hexahedron';
[filename, pathname] = uigetfile( ...
       {'*.dat';'*.txt';'*.*'}, ...
        'Pick a file');
name=[pathname,filename];
%     [type_2d,element,nod,dir,nxe,nye,nip
fid=fopen(name,'rt');
nod=sscanf(fgetl(fid),'%d');
tmp=sscanf(fgetl(fid),'%d %d %d %d %d');
nxe=tmp(1);nye=tmp(2);nze=tmp(3);nip=tmp(4);np_types=tmp(5);
prop=zeros(2,np_types);
for i=1:np_types
    tmp=sscanf(fgetl(fid),'%f %f');
    prop(:,i)=tmp;
end
[nels,nn] = mesh_size( element,nod );
etype=ones(nels,1);
if np_types>1
   for i=1:nels
       etype(i)=sscanf(fgetl(fid),'%d');
   end
end
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
str='z_coords=sscanf(fgetl(fid),''';
for i=1:nze
    str=[str,'%f '];
end
str=[str,''');'];
eval(str);
% x_coords=sscanf(fgetl(fid),'%f %f %f');
% y_coords=sscanf(fgetl(fid),'%f %f %f');
nr=sscanf(fgetl(fid),'%d');
nf=ones(nodof,nn);
for i=1:nr
    str=fgetl(fid);
    tmp=sscanf(str,'%d %d %d %d');
    nf(:,tmp(1))=[tmp(2),tmp(3),tmp(4)];
end
loaded_nodes=sscanf(fgetl(fid),'%d');
if loaded_nodes>0
    loadstmp=zeros(loaded_nodes,4);
    for i=1:loaded_nodes
        tmp=sscanf(fgetl(fid),'%d %f %f %f');
        loadstmp(i,:)=tmp;
    end
end
fixed_freedoms=sscanf(fgetl(fid),'%d');
if fixed_freedoms~=0
    fixedtmp=zeros(fixed_freedoms,3);
    for i=1:fixed_freedoms
        tmp=sscanf(fgetl(fid),'%d %d %f');
        fixedtmp(i,:)=tmp;
    end
end
fclose(fid);

% % nst=6; %The row number of B Matrix
% if strcmp(type_2d,'axisymmetric')
%     nst=4;
% end
[ nels,nn ] = mesh_size( element,nod );
% nf2=ones(nodof,nn);
% for i=1:nr
%     nf2(:,nf(i))=[nf(i,2);nf(i,3)];
% end
% nf=nf2;
[ nf ] = formnf( nf );
neq=max(nf(:));kdiag=zeros(neq,1);
ndof=nod*nodof;
g_num=zeros(nod,nels);
g_coord=zeros(ndim,nn);
g_g=zeros(ndof,nels);
g=zeros(ndof,1); %每个单元的自由度序号

for iel=1:nels
    [ coord,num] = hexahedron_xz(iel,x_coords,y_coords,z_coords,nod);
%     [ num,coord ] = geom_rect(iel,element,x_coords,y_coords,z_coords,nod ,dir);
    [ g ] = num_to_g(num,nf);
    g_num(:,iel)=num;
    g_coord(:,num)=coord';
    g_g(:,iel)=g ;
    kdiag = fkdiag(kdiag, g);
end
% % Mesh wirte into ps 
% % 暂时不要，需要的话后面可以添加
for i=2:neq
    kdiag(i)=kdiag(i)+kdiag(i-1);
end
fid=fopen('result.dat','wt');
str=['There are ',num2str(neq) ' ,equations and the skyline storage is ',num2str(kdiag(neq))];
fprintf(fid,'%s \n',str);
disp(['There are ',num2str(neq) ' ,equations and the skyline storage is ',num2str(kdiag(neq))]);

%% element stiffness integration and assembly
[points, weights]= sample(element,nip );
ih=nst;%ih=3 (plane strain)
%ih=4 (axisymmetry or plane strain elastoplasticity) 
% ih=6(three dimensions)

gc=1;
kv=zeros(kdiag(neq),1);
for iel=1:nels
    dee=deemat(prop(1,etype(iel)),prop(2,etype(iel)),ih);
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
%         deriv=jac/der;
        jac=inv(jac);
        deriv=jac*der;
        bee=beemat( deriv,ih );
%         if strcmp(type_2d,'axisymmetric')
%             gc=fun*coord;
% %             t=1:2:ndof-1;
%             bee(4,1:2:ndof-1)=fun/gc(1);
%         end
        km=km+bee'*dee*bee*Det*weights(i)*gc(1);
    end
    kv=fsparv(kv,km,g,kdiag);
end
% Form the load vector
loads=zeros(neq,1);
if loaded_nodes>0
for i=1:size(loadstmp,1)
    index=nf(:,loadstmp(i,1));
    for j=1:3
        if index(j)~=0
            loads(index(j))=loadstmp(i,j+1);
        end
    end
end
end
if fixed_freedoms ~=0
    node=fixedtmp(:,1);
    sense=fixedtmp(:,2);
    value=fixedtmp(:,3);
    for i=1:fixed_freedoms
        no=nf(sense(i),node(i));
        kv(kdiag(no))=kv(kdiag(no))+penalty;
        loads(no)=kv(kdiag(no))*value(i);
    end    
end
% !-----------------------equation solution---------------------------------
kv=sparin(kv,kdiag);
loads = spabac(kv,loads,kdiag);
disp=zeros(nn,nodof);
fprintf(fid,'Node   x-disp  y-disp \n');
for i=1:nn
    fprintf(fid,' %d ',i);
    for j=1:nodof
        if nf(j,i)==0
            disp(i,j)=0;
        else
            disp(i,j)=loads(nf(j,i));
        end
        fprintf(fid,'%10.8f ',disp(i,j));
    end
    fprintf(fid,'\n');
end

% !-----------------------recover stresses at nip integrating points--------
nip=1; %积分点数量改为1 单元应力
[points, weights]= sample(element,nip );
fprintf(fid,'The integration point stresses are:\n');
fprintf(fid,'Element x-coord     y-coord     z-coord    sig_x   sig_ysig_y   tau_xy   tau_yz   tau_zx\n');
for iel=1:nels
     dee=deemat(prop(1,etype(iel)),prop(2,etype(iel)),ih);
     num=g_num(:,iel);
     coord=g_coord(:,num)';
     g=g_g(:,iel);
     eld=zeros(length(g),1);
     for i=1:length(g)
         if g(i)==0
             eld(i)=0;
         else
             eld(i)=loads(g(i));
         end
     end
     for i=1:nip
         fun=shape_fun(points,i,ndim,nod);
         der=shape_der(points,i,ndim,nod);
         gc=fun'*coord; %高斯点坐标
         jac=der*coord;
         jac=inv(jac);
         deriv=jac*der;
         bee=beemat( deriv,ih );
        strain=bee*eld;
        sigma=dee*(bee*eld);
        fprintf(fid,'%d %10.8f %10.8f  %10.8f %10.8f %10.8f %10.8f %10.8f %10.8f %10.8f\n',iel,gc(1),gc(2),gc(3),...
            sigma(1),sigma(2),sigma(3),sigma(4),sigma(5),sigma(6));
     end
end
fclose(fid);
