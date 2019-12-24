% 2D Viscoplastic of slope
clear;clc;
addpath ../geom
addpath ../main
zero=0;two=2; dt=1.0e15;d3=3.0;d4=4.0;d6=6;one=1.0;%nprops=6;
penalty=1e20;start_dt=1.e15;nip=4;nst=4;
element='quadrilateral';nod=8;ndim=2;nodof=2;ndof=nodof*nod;
[filename, pathname] = uigetfile( ...
    {'*.dat';'*.txt';'*.*'}, ...
    'Pick a file');
name=[pathname,filename];
tlabel=filename(6:length(filename)-4);
%     [type_2d,element,nod,dir,nxe,nye,nip
fid=fopen(name,'rt');
tmp=sscanf(fgetl(fid),'%f %f %f %f %f'); %w1 s1 w2 s2 h1 h2
w1=tmp(1);s1=tmp(2);w2=tmp(3);h1=tmp(4);h2=tmp(5);
tmp=sscanf(fgetl(fid),'%f %f %f %f'); %nx1 nx2 ny1 ny2
nx1=tmp(1);nx2=tmp(2);
ny1=tmp(3);ny2=tmp(4);
np_types=sscanf(fgetl(fid),'%d');
prop=sscanf(fgetl(fid),'%f %f %f %f %f'); %phi c psi gamma e v
tmp=sscanf(fgetl(fid),'%f %f ');
tol=tmp(1);limit=tmp(2);
nsrf=sscanf(fgetl(fid),'%d');
srf=sscanf(fgetl(fid),'%f %f %f %f %f %f');
fclose(fid);
%
nye=ny1+ny2;
nels=nx1*nye+ny2*nx2;
nn=(3*nye+2)*nx1+2*nye+1+(3*ny2+2)*nx2;
nf = emb_2d_bc(nx1,nx2,ny1,ny2);
neq=max(nf(:));

% loop to get elements and nodes
etype=ones(nels,1); %
kdiag=zeros(neq,1);
g_num=zeros(nod,nels);
g_coord=zeros(ndim,nn);
g_g=zeros(ndof,nels);devp=zeros(4,1);
g=zeros(ndof,1); % The global dof number
for iel=1:nels
    [coord,num]=emb_2d_geom(iel,nx1,nx2,ny1,ny2,w1,s1,w2,h1,h2);
    g_num(:,iel)=num;
%     plot(coord(:,1),coord(:,2),'-*');hold on
    [ g ] = num_to_g(num,nf);
    g_coord(:,num)=coord';
    g_g(:,iel)=g ;
    kdiag = fkdiag(kdiag, g);
end
for i=2:neq
    kdiag(i)=kdiag(i)+kdiag(i-1);
end
disp(['There are ',num2str(neq) ' equations and the skyline storage is ',num2str(kdiag(neq))]);
[points, weights]= sample(element,nip );
kv=zeros(kdiag(neq),1);
gravlo=zeros(neq,1); % gravity
% % !-----------------------element stiffness integration and assembly--------
for iel=1:nels
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
end
% !-----------------------factorise equations-------------------------------
kv=sparin(kv,kdiag);
disp([' srf       max disp    iters']);
for iy=1:nsrf
    dt=start_dt;
    for i=1:np_types
        phi=prop(1,i);
        tnph=tan(phi*pi/180);
        phif=atan(tnph/srf(iy));
        snph=sin(phif) ;
        e=prop(5,i) ;
        v=prop(6,i);
        ddt=4*(1+v)*(1-2*v)/(e*(1-2*v+snph^2));
        if ddt<dt
            dt=ddt;
        end
    end
    iters=0;
    bdylds=zeros(neq,1);
    evpt=zeros(nst,nip,nels);
    oldis=zeros(neq,1);
    %-----------------------plastic iteration loop----------------------------
    while 1
        fmax=0;
        iters=iters+1 ;
        loads=gravlo+bdylds ;
        loads = spabac(kv,loads,kdiag);
        tensor=zeros(nst,nip,nels);
        if (iy==1 && iters==1)
            elastic=loads;
        end
        % !-----------------------check plastic convergence-------------------------
        [converged,oldis,res] = checon(loads,oldis,tol);
        if iters==1
            converged=0;
        end
        if converged==1 || iters==limit
            bdylds=zeros(neq,1);
        end
        % !-----------------------go round the Gauss Points -----------------------
        for iel=1:nels %elements_3
            bload=zeros(ndof,1);
            phi=prop(1,etype(iel));
            tnph=tan(phi*pi/180);
            phif=atan(tnph/srf(iy))*180/pi ;
            psi=prop(3,etype(iel));
            tnps=tan(psi*pi/180) ;
            psif=atan(tnps/srf(iy))*180/pi;
            cf=prop(2,etype(iel))/srf(iy) ;
            e=prop(5,etype(iel)) ;
            v=prop(6,etype(iel)) ;
            dee = deemat(e,v,nst);
            
%             c=prop(2,etype(iel));
%             psi=prop(3,etype(iel)) ;
%             dee=deemat(prop(5,etype(iel)),prop(6,etype(iel)),nst);
            
            rload=zeros(ndof,1);% 
            num=g_num(:,iel);
            coord=g_coord(:,num)';
            g=g_g(:,iel);
            eld=zeros(length(g),1);
            for ii=1:length(g) %eld=loads(g) 
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
                f=mocouf(phif,cf,sigm,dsbar,lode_theta);
                if converged==1 && iters==limit
                    devp=stress;
                else
                    if f>=zero
                        [dq1,dq2,dq3]=mocouq(psif,dsbar,lode_theta);
                        [m1,m2,m3]= formm(stress);
                        flow=f*(m1*dq1+m2*dq2+m3*dq3);
                        erate=flow*stress;
                        evp=erate*dt;
                        evpt(:,i,iel)=evpt(:,i,iel)+evp;
                        devp=dee*evp;
                    end
                end
                if f>=0 %|| (converged==1 || iters==limit)
                    eload=bee'*devp;
                    bload=bload+eload*Det*weights(i);
                end
            end
            % !-----------------------compute the total bodyloads vector----------------
            for jj=1:length(g)
                if g(jj)~=0
                    bdylds(g(jj))=bdylds(g(jj))+bload(jj);
            % react(g(jj))=react(g(jj))+rload(jj);
                end
            end
        end %elements_3
        if converged==1 || iters==limit
            break;
        end
    end
    disp([num2str(srf(iy)),'     ',num2str(max(abs(loads))),'         ',num2str(iters)]);
    if iters==limit
        break;
    end
end

% % !-----------------------recover stresses at nip integrating points--------
% nip=4; %���ֵ�������Ϊ1 ��ԪӦ��
% [points, weights]= sample(element,nip );
% 
% % fprintf(fid,'The integration point stresses are:\n');
% % if ndim==2
% %     fprintf(fid,'Element x-coord     y-coord    sig_x   sig_y  tau_xy\n');
% % else
% %     fprintf(fid,'Element x-coord     y-coord     z-coord    sig_x   sig_ysig_y   tau_xy   tau_yz   tau_zx\n');
% % end
% if ndim==2
%     ih=3;
% else
%     ih=6;
% end
% ele_stress=zeros(nels,ih);
% ele_strain=zeros(nels,ih);
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
%      avg_stress=0;avg_strain=0;
%      for i=1:nip
%          fun=shape_fun(points,i,ndim,nod);
%          der=shape_der(points,i,ndim,nod);
%          gc=fun'*coord; %��˹������
%          jac=der*coord;
%          deriv=jac\der;
% %          jac=inv(jac);
% %          deriv=jac*der;
%          bee=beemat( deriv,ih );
%         strain=bee*eld;
%         sigma=dee*(bee*eld);
%         avg_stress=avg_stress+sigma;
%         avg_strain=avg_strain+strain;
% %         if ndim==2
% %             fprintf(fid,'%d %10.8g %10.8g %10.8g %10.8g %10.8g \n',iel,gc(1),gc(2),sigma(1),sigma(2),sigma(3));
% %         else
% %             fprintf(fid,'%d %10.8g %10.8g  %10.8g %10.8g %10.8g %10.8g %10.8g %10.8g %10.8g\n',iel,gc(1),gc(2),gc(3),...
% %                 sigma(1),sigma(2),sigma(3),sigma(4),sigma(5),sigma(6));
% %         end
%      end
%      ele_stress(iel,:)=avg_stress/nip;
%      ele_strain(iel,:)=avg_strain/nip;
% end
% nod_stress=zeros(nn,3);
% nod_strain=zeros(nn,3);
% for k = 1:nn
%     sigx = 0. ;sigy = 0.; tau = 0.;
%     ne = 0;  
%     for iel = 1:nels
%         num=g_num(:,iel);
%         for i=1:length(num)
%             if num(i) == k
%                 ne=ne+1;
%                 nod_stress(k,:)=nod_stress(k,:)+ele_stress(iel,:);
%                 nod_strain(k,:)=nod_strain(k,:)+ele_strain(iel,:);
%             end
%         end
%     end
%     nod_stress(k,:)=nod_stress(k,:)/ne;
%     nod_strain(k,:)=nod_strain(k,:)/ne;
% end
% % fclose(fid);
% mesh_ensi(filename,ndim,nn,nod,element,nels,g_coord,g_num,etype,nf,gravlo,loads,nod_strain,nod_stress);

% % !-----------------------recover stresses at nip integrating points around gauss point--------
% nip=4; %���ֵ�������Ϊ1 ��ԪӦ��
% [points, weights]= sample(element,nip );
% if ndim==2
%     ih=3;
% else
%     ih=6;
% end
% nod_stress=zeros(3,nn);
% nod_strain=zeros(3,nn);
% num_repeat=zeros(1,nn);
% fn=filename(1:length(filename)-4);fcase=[fn,'gauss.dat'];
% fid=fopen(fcase,'wt');
% for iel=1:nels
%      dee=deemat(prop(5,etype(iel)),prop(6,etype(iel)),ih);
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
%      ele_stress=zeros(ih,nip);ele_strain=zeros(ih,nip);
%      for i=1:nip
%          fun=shape_fun(points,i,ndim,nod);
%          der=shape_der(points,i,ndim,nod);
%          gc=fun'*coord; %��˹������
%          jac=der*coord;
%          deriv=jac\der;
%          bee=beemat( deriv,ih );
%          strain=bee*eld;
%          bee2=beemat2( deriv,ih );
%          strain2=bee2*eld;
%          fprintf(fid,'%f  %f  %f  %f\n',strain2(1),strain2(2),strain2(3),strain2(4));
%          sigma=dee*(bee*eld);
%          ele_stress(:,i)=sigma;
%          ele_strain(:,i)=strain;
%          if iel==400 && i==4
%              sigmatmp=sigma;
%              straintmp=strain;
%          end
%      end
%      nod_stress(:,num(1))=nod_stress(:,num(1))+ele_stress(:,3);
%      nod_stress(:,num(2))=nod_stress(:,num(2))+(ele_stress(:,3)+ele_stress(:,1))/2;
%      nod_stress(:,num(3))=nod_stress(:,num(3))+ele_stress(:,1);
%      nod_stress(:,num(4))=nod_stress(:,num(4))+(ele_stress(:,2)+ele_stress(:,1))/2;
%      nod_stress(:,num(5))=nod_stress(:,num(5))+ele_stress(:,2);
%      nod_stress(:,num(6))=nod_stress(:,num(6))+(ele_stress(:,2)+ele_stress(:,4))/2;
%      nod_stress(:,num(7))=nod_stress(:,num(7))+ele_stress(:,4);
%      nod_stress(:,num(8))=nod_stress(:,num(8))+(ele_stress(:,4)+ele_stress(:,3))/2;
%      
%      nod_strain(:,num(1))=nod_strain(:,num(1))+ele_strain(:,3);
%      nod_strain(:,num(2))=nod_strain(:,num(2))+(ele_strain(:,3)+ele_strain(:,1))/2;
%      nod_strain(:,num(3))=nod_strain(:,num(3))+ele_strain(:,1);
%      nod_strain(:,num(4))=nod_strain(:,num(4))+(ele_strain(:,2)+ele_strain(:,1))/2;
%      nod_strain(:,num(5))=nod_strain(:,num(5))+ele_strain(:,2);
%      nod_strain(:,num(6))=nod_strain(:,num(6))+(ele_strain(:,2)+ele_strain(:,4))/2;
%      nod_strain(:,num(7))=nod_strain(:,num(7))+ele_strain(:,4);
%      nod_strain(:,num(8))=nod_strain(:,num(8))+(ele_strain(:,4)+ele_strain(:,3))/2;
%      for i=1:length(num)    
%       num_repeat(num(i))=num_repeat(num(i))+1;
%      end
% end
% fclose(fid);
% nod_stress=nod_stress./num_repeat;
% nod_strain=nod_strain./num_repeat;
% % mesh_ensi(filename,ndim,nn,nod,element,nels,g_coord,g_num,etype,nf,gravlo,loads,nod_strain',nod_stress');
% 
% % %
% % nk=0;
% % for iel=1:nels
% % num=g_num(:,iel);
% % coord=g_coord(:,num)';
% % plot(coord(:,1),coord(:,2),'b-');hold on
% % center=sum(coord)/8;
% % text(center(1)-0.15,center(2),num2str(iel),'fontsize',8);
% % gpos=zeros(nip,2);
% % for i=1:nip
% %     nk=nk+1;
% %     fun=shape_fun(points,i,ndim,nod);
% %     for j=1:2
% %         gpos(i,j)=coord(:,j)'*fun;
% %     end
% % %     plot(gpos(i,1),gpos(i,2),'b*');
% %     text(gpos(i,1),gpos(i,2),num2str(nk),'fontsize',8);
% % end
% % plot(gpos(:,1),gpos(:,2),'b*');
% % end
% 
% 
% 
% % !-----------------------recover debond at nip integrating points around gauss point--------
% nip=4;% Gauss Points Number
% nod_debond=zeros(1,nn);
% nod_rot=zeros(1,nn);
% nod_fabric=zeros(3,nn);
% num_repeat=zeros(1,nn);
% fcase=['slope',tlabel,'.ensi.','debond'];
% dbond=load(['debond',tlabel,'.dat']);
% drot=load(['rot',tlabel,'.dat']);
% fabric=load(['fabric',tlabel,'.dat']);
% detf=[ones(nels*nip,1),zeros(nels*nip,1),zeros(nels*nip,1),ones(nels*nip,1),ones(nels*nip,1),zeros(nels*nip,1),zeros(nels*nip,1),ones(nels*nip,1),ones(nels*nip,1),zeros(nels*nip,1),zeros(nels*nip,1),ones(nels*nip,1)]/2;
% tfabric=(fabric-detf)*4;
% tfab=tfabric.*tfabric;
% dfabric=[sqrt(tfab(:,1)+tfab(:,2)+tfab(:,3)+tfab(:,4)),sqrt(tfab(:,5)+tfab(:,6)+tfab(:,7)+tfab(:,8)),sqrt(tfab(:,9)+tfab(:,10)+tfab(:,11)+tfab(:,12))]/sqrt(2);
% fid=fopen(fcase,'wt');
% k=0;
% for iel=1:nels
%     ele_debond=zeros(1,nip);
%     ele_rot=zeros(1,nip);
%     ele_fabric=zeros(3,nip);
%     num=g_num(:,iel);
%     g=g_g(:,iel);
%      for i=1:nip
%          k=k+1;
%          ele_debond(:,i)=dbond(k);
%          ele_rot(:,i)=drot(k);
%          ele_fabric(:,i)=dfabric(k,:)';
%          
%      end
%      nod_debond(:,num(1))=nod_debond(:,num(1))+ele_debond(:,3);
%      nod_debond(:,num(2))=nod_debond(:,num(2))+(ele_debond(:,3)+ele_debond(:,1))/2;
%      nod_debond(:,num(3))=nod_debond(:,num(3))+ele_debond(:,1);
%      nod_debond(:,num(4))=nod_debond(:,num(4))+(ele_debond(:,2)+ele_debond(:,1))/2;
%      nod_debond(:,num(5))=nod_debond(:,num(5))+ele_debond(:,2);
%      nod_debond(:,num(6))=nod_debond(:,num(6))+(ele_debond(:,2)+ele_debond(:,4))/2;
%      nod_debond(:,num(7))=nod_debond(:,num(7))+ele_debond(:,4);
%      nod_debond(:,num(8))=nod_debond(:,num(8))+(ele_debond(:,4)+ele_debond(:,3))/2;
%      
%      nod_rot(:,num(1))=nod_rot(:,num(1))+ele_rot(:,3);
%      nod_rot(:,num(2))=nod_rot(:,num(2))+(ele_rot(:,3)+ele_rot(:,1))/2;
%      nod_rot(:,num(3))=nod_rot(:,num(3))+ele_rot(:,1);
%      nod_rot(:,num(4))=nod_rot(:,num(4))+(ele_rot(:,2)+ele_rot(:,1))/2;
%      nod_rot(:,num(5))=nod_rot(:,num(5))+ele_rot(:,2);
%      nod_rot(:,num(6))=nod_rot(:,num(6))+(ele_rot(:,2)+ele_rot(:,4))/2;
%      nod_rot(:,num(7))=nod_rot(:,num(7))+ele_rot(:,4);
%      nod_rot(:,num(8))=nod_rot(:,num(8))+(ele_rot(:,4)+ele_rot(:,3))/2;
%      
%      nod_fabric(:,num(1))=nod_fabric(:,num(1))+ele_fabric(:,3);
%      nod_fabric(:,num(2))=nod_fabric(:,num(2))+(ele_fabric(:,3)+ele_fabric(:,1))/2;
%      nod_fabric(:,num(3))=nod_fabric(:,num(3))+ele_fabric(:,1);
%      nod_fabric(:,num(4))=nod_fabric(:,num(4))+(ele_fabric(:,2)+ele_fabric(:,1))/2;
%      nod_fabric(:,num(5))=nod_fabric(:,num(5))+ele_fabric(:,2);
%      nod_fabric(:,num(6))=nod_fabric(:,num(6))+(ele_fabric(:,2)+ele_fabric(:,4))/2;
%      nod_fabric(:,num(7))=nod_fabric(:,num(7))+ele_fabric(:,4);
%      nod_fabric(:,num(8))=nod_fabric(:,num(8))+(ele_fabric(:,4)+ele_fabric(:,3))/2;
%      
%      for i=1:length(num)    
%       num_repeat(num(i))=num_repeat(num(i))+1;
%      end
% end
% fclose(fid);
% nod_debond=nod_debond./num_repeat;
% nod_rot=rad2deg(nod_rot./num_repeat);
% nod_fabric=nod_fabric./num_repeat;
% mesh_ensi(filename,ndim,nn,nod,element,nels,g_coord,g_num,etype,nf,gravlo,loads,nod_strain',nod_stress',nod_debond,nod_rot,nod_fabric);
% 
% 
% fid=fopen(['point_out',tlabel,'.dat'],'wt');
% ielt=[400];% 218
% for kk=1
%     iel=ielt(kk);
%     g=g_g(:,iel);
%     eld=zeros(length(g),1);
%     for ii=1:length(g) %eld=loads(g)
%         if g(ii)~=0
%             eld(ii)=loads(g(ii));
%         else
%             eld(ii)=0;
%         end
%     end
%     i=4;
%     fun=shape_fun(points,i,ndim,nod);
%     dpx=fun'*eld([1:2:16]);
%     dpy=fun'*eld([1:2:16]);
%     dp=sqrt(dpx^2+dpy^2);
%     fprintf(fid,'%f\n',dp);
% end
% fcase=[fn,'gauss.dat'];
% tgauss=load(fcase);
% straint=tgauss(1600,:);
% sm=(straint(1)+straint(2))/2;
% shearEt=(2*((straint(1)-sm)^2+(straint(2)-sm)^2+straint(3)^2+straint(4)^2)).^0.5;
% fprintf(fid,'%f\n',shearEt);
% sm=(sigmatmp(1)+sigmatmp(2))/2;
% shearSt=(0.5*((sigmatmp(1)-sm)^2+(sigmatmp(2)-sm)^2+sigmatmp(3)^2+sigmatmp(3)^2)).^0.5;
% fprintf(fid,'%f\n',shearSt);
% fprintf(fid,'%f\n',dbond(1600));
% fprintf(fid,'%f\n',dfabric(1600,1));
% fprintf(fid,'%f\n',rad2deg(drot(1600,1)));
% fclose(fid);
