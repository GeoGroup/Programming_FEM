clear;clc;
addpath ../library/geom
addpath ../library/main
penalty=1e20;
ndim=2;nodof=2;
% global nxe nye nze;
[filename, pathname] = uigetfile( ...
       {'*.dat';'*.txt';'*.*'}, ...
        'Pick a file');
name=[pathname,filename];
%     [type_2d,element,nod,dir,nxe,nye,nip
fid=fopen(name,'rt');
type_2d=fgetl(fid);type_2d([1,length(type_2d)])=[];
element=fgetl(fid);element([1,length(element)])=[];
nod=sscanf(fgetl(fid),'%d');
dir=fgetl(fid);dir([1,length(dir)])=[];
tmp=sscanf(fgetl(fid),'%d %d %d %d');
nxe=tmp(1);
nye=tmp(2);
nip=tmp(3);
np_types=tmp(4);
prop=zeros(4,np_types);
for i=1:np_types
    tprop=sscanf(fgetl(fid),'%f %f');
    prop(:,i)=tprop';
end

if np_types>1
    str=fgetl(fid);
    etype=str2num(cell2mat(split(str,' ')));
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
% x_coords=sscanf(fgetl(fid),'%f %f %f');
% y_coords=sscanf(fgetl(fid),'%f %f %f');

[ nels,nn ] = mesh_size(element,nod,nxe,nye );
if np_types==1
    etype=ones(nels,1);    
end

nr=sscanf(fgetl(fid),'%d');
nf=zeros(nr,3);
for i=1:nr
    tmp=sscanf(fgetl(fid),'%d %d %d');
    nf(i,:)=tmp;
end

dtemp=zeros(nn,1);
id=0;
for i=1:41
    str=fgetl(fid);
    tmp=split(str,' ');
    tmp(cellfun(@isempty,tmp))=[];
    tmp2=zeros(length(tmp),1);
    for j=1:length(tmp)
        tmp2(j)=str2double(tmp{j});
    end        
    dtemp(id+1:id+length(tmp))=tmp2;
    id=id+length(tmp2);
end

nspr=sscanf(fgetl(fid),'%d');
spno=zeros(nspr,1);
spse=zeros(nspr,1);
spva=zeros(nspr,1);
for i=1:nspr
    tmp=sscanf(fgetl(fid),'%d %d %f');
    spno(i)=tmp(1);
    spse(i)=tmp(2);
    spva(i)=tmp(3);
end

loaded_nodes=sscanf(fgetl(fid),'%d');
if loaded_nodes>0
    loadstmp=zeros(loaded_nodes,3);
    for i=1:loaded_nodes
        tmp=sscanf(fgetl(fid),'%d %f %f');
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

nst=3; %The row number of B Matrix
if strcmp(type_2d,'axisymmetric')
    nst=4;
end

nf2=ones(nodof,nn);
for i=1:nr
    nf2(:,nf(i))=[nf(i,2);nf(i,3)];
end
nf=nf2;
[ nf ] = formnf( nf );
neq=max(nf(:));kdiag=zeros(neq,1);
ndof=nod*nodof;
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



% Mesh wirte into ps 
% 暂时不要，需要的话后面可以添加
for i=2:neq
    kdiag(i)=kdiag(i)+kdiag(i-1);
end
fid=fopen('result.dat','wt');
str=['There are ',num2str(neq) ' equations and the skyline storage is ',num2str(kdiag(neq))];
fprintf(fid,'%s \n',str);
disp(['There are ',num2str(neq) ' equations and the skyline storage is ',num2str(kdiag(neq))]);

%% element stiffness integration and assembly
[points, weights]= sample(element,nip );
ih=nst;%ih=3 (plane strain)
%ih=4 (axisymmetry or plane strain elastoplasticity) 
% ih=6(three dimensions)
% etype=ones(nels,1);
gc=1;
tload=zeros(neq,1);
kv=zeros(kdiag(neq),1);
teps=zeros(nst,1);
for iel=1:nels
    dee=deemat(prop(1,etype(iel)),prop(2,etype(iel)),ih);
    num=g_num(:,iel);
    coord=g_coord(:,num)';
    g=g_g(:,iel);
    dtel=dtemp(num);etl=zeros(ndof,1);
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
        if strcmp(type_2d,'axisymmetric')
            gc=fun'*coord;
%             t=1:2:ndof-1;
            bee(4,1:2:ndof-1)=fun'/gc(1);
        end
        km=km+bee'*dee*bee*Det*weights(i)*gc(1);
        gtemp=fun'*dtel;
        teps(1:2)=gtemp*[prop(3,etype(iel)),prop(4,etype(iel))];
        etl=etl+bee'*dee*teps*Det*weights(i)*gc(1);
    end
    for i=1:length(g)
        if g(i)~=0
            tload(g(i))=tload(g(i))+etl(i);
        end
    end
    kv=fsparv(kv,km,g,kdiag);
end
% 
if nspr>0
    sno=zeros(size(nspr));
    for i=1:nspr
        sno(i)=nf(spse(i),spno(i));
    end
    kv(kdiag(sno))=kv(kdiag(sno))+spva;
end

% Form the load vector
loads=zeros(neq,1);
if loaded_nodes>0
for i=1:size(loadstmp,1)
    index=nf(:,loadstmp(i,1));
    for j=1:2
        if index(j)~=0
            loads(index(j))=loadstmp(i,j+1);
        end
    end
end
end
loads=loads+tload;
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

% % !-----------------------recover stresses at nip integrating points--------
% nip=1; %积分点数量改为1 单元应力
% [points, weights]= sample(element,nip );
% fprintf(fid,'The integration point stresses are:\n');
% fprintf(fid,'Element x-coord     y-coord    sig_x   sig_y   tau_xy\n');
% for iel=1:nels
%      dee=deemat(prop(1,etype(iel)),prop(2,etype(iel)),ih);
%      num=g_num(:,iel);
%      coord=g_coord(:,num)';
%      g=g_g(:,iel);
%      eld=zeros(length(g),1);
%      for i=1:length(g)
%          if g(i)==0
%              eld(i)=0;
%          else
%              eld(i)=loads(g(i));
%          end
%      end
%      for i=1:nip
%          fun=shape_fun(points,i,ndim,nod);
%          der=shape_der(points,i,ndim,nod);
%          gc=fun'*coord; %高斯点坐标
%          jac=der*coord;
%          jac=inv(jac);
%          deriv=jac*der;
%          bee=beemat( deriv,ih );
%         if strcmp(type_2d,'axisymmetric')
%             gc=fun'*coord;
% %             t=1:2:ndof-1;
%             bee(4,1:2:ndof-1)=fun/gc(1);
%         end
%         strain=bee*eld;
%         sigma=dee*(bee*eld);
%         fprintf(fid,'%d %10.8f %10.8f %10.8f %10.8f %10.8f \n',iel,gc(1),gc(2),sigma(1),sigma(2),sigma(3));
%      end
% end
% fclose(fid);
