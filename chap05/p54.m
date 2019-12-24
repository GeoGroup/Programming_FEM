clear;clc;
addpath ../geom
addpath ../main
penalty=1e20;
[filename, pathname] = uigetfile( ...
       {'*.dat';'*.txt';'*.*'}, ...
        'Pick a file');
name=[pathname,filename];
%     [type_2d,element,nod,dir,nxe,nye,nip
fid=fopen(name,'rt');
element=fgetl(fid);element([1,length(element)])=[];
tmp=sscanf(fgetl(fid),'%d %d %d %d %d %d %d %d');
nod=tmp(1);nels=tmp(2);
nn=tmp(3);nip=tmp(4);
nodof=tmp(5);nst=tmp(6);
ndim=tmp(7);np_types=tmp(8);
ndof=nod*nodof;
prop=zeros(3,np_types);
for i=1:np_types
    str=fgetl(fid);
    prop(:,i)=sscanf(str,'%f %f %f');
end
etype=ones(nels,1);
if np_types>1
    etype=sscanf(fgetl(fid),'%d');
end
g_coord=zeros(ndim,nn);
for i=1:nn
    if ndim==2
        g_coord(:,i)=sscanf(fgetl(fid),'%f %f');
    else
        g_coord(:,i)=sscanf(fgetl(fid),'%f %f %f');
    end
end
g_num=zeros(nod,nels);
str='g_num(:,i)=sscanf(fgetl(fid),''';
for i=1:nod
    str=[str,'%f '];
end
str=[str,''');'];
for i=1:nels
    eval(str);
end
nf=ones(nodof,nn); %
nr=sscanf(fgetl(fid),'%d');
for i=1:nr
    if ndim==2
        tmp=sscanf(fgetl(fid),'%d %d %d');
    else
        tmp=sscanf(fgetl(fid),'%d %d %d %d');
    end
    nf(:,tmp(1))=tmp(2:length(tmp));
end
nf=formnf(nf);
neq=max(nf(:));
kdiag=zeros(neq,1);%
g_g=zeros(ndof,nels); %
%% ! -----------------------loop the elements to find global arrays sizes-----
for iel=1:nels
    num=g_num(:,iel);
    [ g ] = num_to_g(num,nf);
    g_num(:,iel)=num;
    g_g(:,iel)=g ;
    kdiag = fkdiag(kdiag, g);
end
for i=2:neq
    kdiag(i)=kdiag(i)+kdiag(i-1);
end
kv=zeros(kdiag(neq));
disp(['There are ',num2str(neq) ' ,equations and the skyline storage is ',num2str(kdiag(neq))]);
%% !-----------------------element stiffness integration and assembly--------
[points, weights]= sample(element,nip );
if ndim==2
    ih=3;
else
    ih=6;
end
kv=zeros(kdiag(neq),1);    
gravlo=zeros(neq,1);
for iel=1:nels
    eld=zeros(ndof,1);
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
        deriv=jac\der;
%         jac=inv(jac);
%         deriv=jac*der;
        bee=beemat( deriv,ih );
        km=km+bee'*dee*bee*Det*weights(i);
        eld(nodof:nodof:ndof)=eld(nodof:nodof:ndof)+fun(:)*Det*weights(i);
    end
    kv=fsparv(kv,km,g,kdiag);
    for i=1:ndof
        if g(i)~=0
            gravlo(g(i))=gravlo(g(i))-eld(i)*prop(3,etype(iel));
        end
    end
end
loaded_nodes=sscanf(fgetl(fid),'%d');
loads=zeros(neq,1);
if loaded_nodes>0
    for i=1:loaded_nodes
        if ndim==2
            tmp=sscanf(fgetl(fid),'%d %f %f');
        else
            tmp=sscanf(fgetl(fid),'%d %f %f %f');
        end
        index=nf(:,tmp(1));
        for j=1:ndim
            if index(j)~=0
                loads(index(j))=tmp(j+1);
            end
        end
    end
end
loads=loads+gravlo;
fixed_freedoms=sscanf(fgetl(fid),'%d');
if fixed_freedoms >0
    for i=1:fixed_freedoms
        tmp=sscanf(fgetl(fid),'%d %d %f');
        node=tmp(1);sense=tmp(2);value=tmp(3);
        no=nf(sense,node);
        kv(kdiag(no))=kv(kdiag(no))+penalty;
        loads(no)=kv(kdiag(no))*value;
    end
end
fclose(fid);
fid=fopen('result.dat','wt');
% !-----------------------equation solution---------------------------------
kv=sparin(kv,kdiag);
loads = spabac(kv,loads,kdiag);
disp=zeros(nn,nodof);loads(abs(loads)<1e-10)=0;
if ndim==2
    fprintf(fid,'Node   x-disp  y-disp \n');
else
    fprintf(fid,'Node   x-disp  y-disp z-disp\n');
end
for i=1:nn
    fprintf(fid,' %d ',i);
    for j=1:nodof
        if nf(j,i)==0
            disp(i,j)=0;
        else
            disp(i,j)=loads(nf(j,i));
        end
        fprintf(fid,'%10.8g ',disp(i,j));
    end
    fprintf(fid,'\n');
end
fclose(fid);
% % !-----------------------recover stresses at nip integrating points--------
% nip=1; %积分点数量改为1 单元应力
% [points, weights]= sample(element,nip );
% 
% fprintf(fid,'The integration point stresses are:\n');
% if ndim==2
%     fprintf(fid,'Element x-coord     y-coord    sig_x   sig_y  tau_xy\n');
% else
%     fprintf(fid,'Element x-coord     y-coord     z-coord    sig_x   sig_ysig_y   tau_xy   tau_yz   tau_zx\n');
% end
% if ndim==2
%     ih=3;
% else
%     ih=6;
% end
% ele_stress=zeros(nels,ih);
% 
% for iel=1:nels
%     
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
%      avg_stress=0;
%      for i=1:nip
%          fun=shape_fun(points,i,ndim,nod);
%          der=shape_der(points,i,ndim,nod);
%          gc=fun'*coord; %高斯点坐标
%          jac=der*coord;
%          deriv=jac\der;
% %          jac=inv(jac);
% %          deriv=jac*der;
%          bee=beemat( deriv,ih );
%         strain=bee*eld;
%         sigma=dee*(bee*eld);
%         avg_stress=avg_stress+sigma;
%         if ndim==2
%             fprintf(fid,'%d %10.8g %10.8g %10.8g %10.8g %10.8g \n',iel,gc(1),gc(2),sigma(1),sigma(2),sigma(3));
%         else
%             fprintf(fid,'%d %10.8g %10.8g  %10.8g %10.8g %10.8g %10.8g %10.8g %10.8g %10.8g\n',iel,gc(1),gc(2),gc(3),...
%                 sigma(1),sigma(2),sigma(3),sigma(4),sigma(5),sigma(6));
%         end
%      end
%      ele_stress(iel,:)=avg_stress/nip;
% end
% nod_stress=zeros(nn,ih);
% for k = 1:nn
%     sigx = 0. ;sigy = 0.; tau = 0.;
%     ne = 0;  
%     for iel = 1:nels
%         num=g_num(:,iel);
%         for i=1:length(num)
%             if num(i) == k
%                 ne=ne+1;
%                 nod_stress(k,:)=nod_stress(k,:)+ele_stress(iel,:);
%             end
%         end
%     end
%     nod_stress(k,:)=nod_stress(k,:)/ne;
% end
% fclose(fid);
