% 2D Viscoplastic 
clear;clc;
addpath ../geom
addpath ../main
zero=0;two=2;dt=1.0e15;d3=3.0;d4=4.0;d6=6;one=1.0;%nprops=6;
penalty=1e20;start_dt=1.e15;pt5=0.5;
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
np_types=tmp(3);
[ nels,nn ] = mesh_size( element,nod );
prop=sscanf(fgetl(fid),'%f %f %f %f %f %f %f');


etype=ones(nels,1);
if np_types>1
    for i=1:nels
        etype(i)=sscanf(fgetl(fid),'%d');
    end
end
% qs=sscanf(fgetl(fid),'%f');
% nf=bc_rect(nxe,nye,'y');
% neq=max(nf(:));

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
    nf(:,tmp(1))=tmp(2:length(tmp));
end
tensor=zeros(nst,nip,nels);
nf=formnf(nf);
neq=max(nf(:));
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

[points, weights]= sample(element,nip );
kv=zeros(kdiag(neq),1);    
% gravlo=zeros(neq,1);
% !-----------------------element stiffness integration and assembly--------
for iel=1:nels % elements_2
    dee=deemat(prop(6,etype(iel)),prop(7,etype(iel)),nst);
    gamma=prop(4,etype(iel));
    k0=prop(5,etype(iel));
    num=g_num(:,iel);
    coord=g_coord(:,num)';
    g=g_g(:,iel);
    km=zeros(ndof,ndof);
%     eld=zeros(ndof,1);
    
    for i=1:nip  % gauss_pts_1
        fun=shape_fun(points,i,ndim,nod);
        der=shape_der(points,i,ndim,nod);
        % ADD Initial Soil pressure
        gc=fun'*coord;
        tensor(2,i,iel)=(gc(2)-y_coords(1))*gamma;
        tensor(1,i,iel)=(gc(2)-y_coords(1))*gamma*k0;
        tensor(4,i,iel)=tensor(1,i,iel);
        jac=der*coord;
        Det=det(jac);
        deriv=jac\der;
        bee=beemat( deriv,nst );
        km=km+bee'*dee*bee*Det*weights(i);
%         eld(nodof:nodof:ndof)=eld(nodof:nodof:ndof)+fun(:)*Det*weights(i);
    end % gauss_pts_1
    kv=fsparv(kv,km,g,kdiag);
end % elements_2

fixed_freedoms=sscanf(fgetl(fid),'%d');
if fixed_freedoms>0
    node=zeros(fixed_freedoms,1);
    sense=zeros(fixed_freedoms,1);
    no=zeros(fixed_freedoms,1);
    storkv=zeros(fixed_freedoms,1);
    for i=1:fixed_freedoms
        tmp=sscanf(fgetl(fid),'%d %d');
        node(i)=tmp(1);sense(i)=tmp(2);
        no(i)=nf(sense(i),node(i));
        kv(kdiag(no(i)))=kv(kdiag(no(i)))+penalty;
        storkv(i)=kv(kdiag(no(i)));
    end
end
kv=sparin(kv,kdiag);

% !-----------------------displacement increment loop-------------------------------
tmp=sscanf(fgetl(fid),'%f %d %d %f');
tol=tmp(1);limit=tmp(2);
incs=tmp(3);presc=tmp(4);
oldis=zeros(neq,1);totd=zeros(neq,1);bdylds=zeros(neq,1);
% devp=zeros(nst,1);
for iy=1:5%incs %disp_incs
%     bdylds=zeros(neq,1);
 iters=0;
react=zeros(neq,1);
   
    evpt=zeros(nst,nip,nels);
%     % !-----------------------plastic iteration loop----------------------------
    while 1 %its
        iters=iters+1;
        loads=bdylds;
%         loads=loads+bdylds;
        for i=1:fixed_freedoms
            loads(no(i))=storkv(i)*presc;
        end

        loads = spabac(kv,loads,kdiag);
        save 'load.txt' loads -ascii
        bdylds=zeros(neq,1);
        % !-----------------------check plastic convergence-------------------------
        [converged,oldis,res] = checon(loads,oldis,tol);  
% %         disp(['res equal to ',num2str(res)]);
        if iters==1
            converged=0;
        end
%         if converged==1 || iters==limit
%             bdylds=zeros(neq,1);
%         end
        % !-----------------------go round the Gauss Points ------------------------
        for iel=1:nels %elements_3
            phi=prop(1,etype(iel));
            c=prop(2,etype(iel));
            psi=prop(3,etype(iel)) ;
            dee=deemat(prop(6,etype(iel)),prop(7,etype(iel)),nst);
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
            for i=1:nip %gauss_pts_2
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
                fnew=mocouf(phi,c,sigm,dsbar,lode_theta);
                elso=zeros(nst,1);
                if fnew>zero
                    stress=tensor(:,i,iel);
                    [sigm,dsbar,lode_theta]=invar(stress);
                    f=mocouf(phi,c,sigm,dsbar,lode_theta);
                    fac=fnew/(fnew-f);
                    stress=(one-fac)*sigma+tensor(:,i,iel);
                    pl=mcdpl(phi,psi,dee,stress);
                    pl=fac*pl;
                    elso=pl*eps;
                    eload=bee*elso;
                    bload=bload+eload*Det*weights(i);
                end
                % !-----------------------update the Gauss Point stresses-------------------
                if converged ==1 || iters==limit
                    tensor(:,i,iel)=tensor(:,i,iel)+stress;
                    rload=rload+bee'*stress*Det*weights(i);
                end
            end %gauss_pts_2
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
        disp([num2str(iters),'   iters   ',num2str(res)]);
    end %its
%     
    totd=totd+loads;
    pr=zero;
    ot=zero;
    pav=zero;
    for i=1:fixed_freedoms
        pr=pr+react(no(i));
        ot=ot+react(no(i))*g_coord(2,node(i));
    end
    for i=1:4
        pav=pav+(y_coords(i)-y_coords(i+1))*(tensor(1,1,i)+tensor(1,3,i))*pt5;
    end
%     pav=pav/(two*nbo2);
    disp([num2str(iy),'   ',num2str(-totd(1)),'  ',num2str(-pr),'  ',num2str(-pav),'   ',num2str(iters)]);
    if iters==limit
        break;
    end
%     
% %     fprintf('%2d %3d %8.5g %3d \n',iy,qinc(iy),loads(nf(2,node(1))),iters);
% % %     disp([num2str(iy),num2str(ptot),num2str(totd(nf(2,node(1)))),num2str(iters)]);
end %disp_incs
