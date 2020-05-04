% % READ 1/Qp AND Qp/Qs MODELS AND CALCULATE 1/Qk
clear;

pardir='~/GoogleDriveMSU/Work/Lau/Qtomo/3D1Dtomo/v2550h80top2_QsQps/sameQps_svd/';
para='2.0_p0.9';
ifhighres=0;
Vratio_err=0.2;

% % % Input files
if ifhighres>0
    Qpdir='~/GoogleDriveMSU/Work/Lau/Qtomo/3D1Dtomo/v25h30top2_Qp/fcp0.5-20MPa0.27/';
    Qpfl=[Qpdir 'Qinv_7.0e-04_1e+06.p'];
    Qpnodfl=fopen([Qpdir 'Q.grid'],'r');
else
    Qpfl=[pardir 'Qinv.p'];
    Qperrfl=[pardir 'Qinv.p_varm'];
end
Qsfl=[pardir 'Qinv_' para '.s'];
QpQsfl=[pardir 'Qinv_' para '.ps'];
QpQserrfl=[pardir 'Qinv_' para '.ps_varm'];
Vnodfl=fopen([pardir 'lau.vel'],'r');
Qnodfl=fopen([pardir 'Q.grid'],'r');

% % % Output files
if ifhighres==1
    Qkfl=[pardir 'Qinv_' para '_high2low.k'];
    QmQkfl=[pardir 'Qinv_' para '_high2low.mk'];
elseif ifhighres==2
    Qkfl=[pardir 'Qinv_' para '_low2high.k'];
    QmQkfl=[pardir 'Qinv_' para '_low2high.mk'];
else
    Qkfl=[pardir 'Qinv_' para '.k'];
    Qkstdfl=[pardir 'Qinv_' para '.k_varm'];
    QmQkfl=[pardir 'Qinv_' para '.mk'];
end

nQs1D=4;
minQ=0.01;

%% IMPORT Qs NODES, low resolution
temp=textscan(Qnodfl,'');
nx=temp{1}(1);ny=temp{2}(1);nz=temp{3}(1);
reflat=temp{1}(2);reflon=temp{2}(2);beta=temp{3}(2);
lat=temp{1}(3:end-1);lon=temp{2}(3:end-1);
Qz=zeros(1,nz);
for i=1:length(Qz)
    Qz(i)=temp{i}(end);
end
Qx=zeros(1,nx);Qy=zeros(1,ny);
for k=1:length(lat)
    if k==1
        i=1;j=1;
        [Qx(i),Qy(j)]=latlon2xy(reflat,reflon,beta,lat(k),lon(k));
    elseif k<=ny
        j=j+1;
        [tmpx,Qy(j)]=latlon2xy(reflat,reflon,beta,lat(k),lon(k));
    elseif rem(k,ny)==1
        i=i+1;
        [Qx(i),tmpy]=latlon2xy(reflat,reflon,beta,lat(k),lon(k));
    end
end
% Qsx=Qx;Qsy=Qy;Qsz=Qz;
[Qsx,Qsy,Qsz]=meshgrid(Qx,Qy,Qz);
Qsnx=nx;Qsny=ny;Qsnz=nz;

% IMPORT Qs MODEL, low resolution
Qtemp=load(Qsfl);
Qs1D=Qtemp(end-nQs1D+1:end);
Qtemp=Qtemp(1:end-nQs1D);
inode=0;
Qinv=zeros(Qsny,Qsnx,Qsnz);
for i=1:Qsnx
    for j=1:Qsny
        for k=1:Qsnz
            inode=inode+1;
            Qinv(j,i,k)=Qtemp(inode);
        end
    end
end
Qsinv=Qinv;
% scatterQmod(Qsnx,Qsny,Qsnz,Qsinv.*1000,2);

% IMPORT Qp/Qs MODEL, low resolution
Qtemp=load(QpQsfl);
Qtemp=Qtemp(1:end-nQs1D);
inode=0;
Qinv=zeros(Qsny,Qsnx,Qsnz);
for i=1:Qsnx
    for j=1:Qsny
        for k=1:Qsnz
            inode=inode+1;
            Qinv(j,i,k)=Qtemp(inode);
        end
    end
end
QpQs=Qinv;
% scatterQmod(Qsnx,Qsny,Qsnz,Qpinv.*1000,2);
% IMPORT Qp/Qs MODEL STD, low resolution
Qtemp=load(QpQserrfl);
Qtemp=Qtemp(1:end-nQs1D,2);
inode=0;
Qinv=zeros(Qsny,Qsnx,Qsnz);
for i=1:Qsnx
    for j=1:Qsny
        for k=1:Qsnz
            inode=inode+1;
            Qinv(j,i,k)=Qtemp(inode);
        end
    end
end
QpQsstd=Qinv;

%% IMPORT Qp NODES, high resolution
if ifhighres>0
    temp=textscan(Qpnodfl,'');
    nx=temp{1}(1);ny=temp{2}(1);nz=temp{3}(1);
    lat=temp{1}(3:end-1);lon=temp{2}(3:end-1);
    Qz=zeros(1,nz);
    Qz(1)=12;
    for i=1:length(Qz)-1
        Qz(i+1)=temp{i}(end);
    end
    Qx=zeros(1,nx);Qy=zeros(1,ny);
    for k=1:length(lat)
        if k==1
            i=1;j=1;
            [Qx(i),Qy(j)]=latlon2xy(reflat,reflon,beta,lat(k),lon(k));
        elseif k<=ny
            j=j+1;
            [tmpx,Qy(j)]=latlon2xy(reflat,reflon,beta,lat(k),lon(k));
        elseif rem(k,ny)==1
            i=i+1;
            [Qx(i),tmpy]=latlon2xy(reflat,reflon,beta,lat(k),lon(k));
        end
    end
    [Qpx,Qpy,Qpz]=meshgrid(Qx,Qy,Qz);
    Qpnx=nx;Qpny=ny;Qpnz=nz;
    nQp1D=8;
else
    Qpnx=Qsnx;Qpny=Qsny;Qpnz=Qsnz;
    nQp1D=4;
end

% IMPORT Qp MODEL
Qtemp=load(Qpfl);
Qtemp=Qtemp(1:end-nQp1D);
inode=0;
Qinv=zeros(Qpny,Qpnx,Qpnz);
for i=1:Qpnx
    for j=1:Qpny
        for k=1:Qpnz
            inode=inode+1;
            Qinv(j,i,k)=Qtemp(inode);
        end
    end
end
Qpinv=Qinv;
% scatterQmod(Qpnx,Qpny,Qpnz,Qpinv.*1000,2);
% IMPORT Qp MODEL STD
Qtemp=load(Qperrfl);
Qtemp=Qtemp(1:end-nQp1D);
inode=0;
Qinv=zeros(Qpny,Qpnx,Qpnz);
for i=1:Qpnx
    for j=1:Qpny
        for k=1:Qpnz
            inode=inode+1;
            Qinv(j,i,k)=Qtemp(inode);
        end
    end
end
Qpinvstd=Qinv;

%% IMPORT V NODES
temp=textscan(Vnodfl,'');
nx=temp{1}(1);ny=temp{2}(1);nz=temp{3}(1);
reflat=temp{1}(2);reflon=temp{2}(2);beta=temp{3}(2);
lat=temp{1}(3:nx*ny+2);lon=temp{2}(3:nx*ny+2);
Qz=zeros(1,nz);
for i=1:length(Qz)
    Qz(i)=temp{i}(nx*ny+3);
end
Qx=zeros(1,nx);Qy=zeros(1,ny);
for k=1:length(lat)
    if k==1
        i=1;j=1;
        [Qx(i),Qy(j)]=latlon2xy(reflat,reflon,beta,lat(k),lon(k));
    elseif k<=ny
        j=j+1;
        [tmpx,Qy(j)]=latlon2xy(reflat,reflon,beta,lat(k),lon(k));
    elseif rem(k,ny)==1
        i=i+1;
        [Qx(i),tmpy]=latlon2xy(reflat,reflon,beta,lat(k),lon(k));
    end
end
% Qpx=Qx;Qpy=Qy;Qpz=Qz;
[Vx,Vy,Vz]=meshgrid(Qx,Qy,Qz);
Vnx=nx;Vny=ny;Vnz=nz;

%% IMPORT Vp MODEL
Vmod=zeros(ny,nx,nz);
iline=nx*ny+3;
for i=1:nx
    for j=1:ny
        iline=iline+1;
        for k=1:nz
            Vmod(j,i,k)=temp{k}(iline);
        end
    end
end
Vpmod=Vmod;
% scatterQmod(nx,ny,nz,Vpmod,2);

%% IMPORT Vs MODEL
Vmod=zeros(ny,nx,nz);
for i=1:nx
    for j=1:ny
        iline=iline+1;
        for k=1:nz
            Vmod(j,i,k)=temp{k}(iline);
        end
    end
end
Vsmod=Vmod;
% scatterQmod(nx,ny,nz,Vsmod,2);

%% Match resolution, downgrade Qp or upgrade Qs
if ifhighres==1  % Downgrade Qp
    % APPLY SMOOTHING ON Qp MODEL
    Qpinvsm=smooth3(Qpinv,'gaussian',ceil(size(Qpinv,1)/size(Qsinv,1))*2+1);
    % INTERPOLATE Qp TO Qs GRID
    Qpinvnew=interp3(Qpx,Qpy,Qpz,Qpinvsm,Qsx,Qsy,Qsz,'spline');
    Qpinv=Qpinvnew;
    imax=Qsnx;jmax=Qsny;kmax=Qsnz;
    Qx=Qsx;Qy=Qsy;Qz=Qsz;
elseif ifhighres==2  % Upgrade Qs
    % INTERPOLATE Qs TO Qp GRID
%     scatterQmod(Qsnx,Qsny,Qsnz,Qsinv,2);caxis([0 0.05]);
    Qsinvnew=interp3(Qsx,Qsy,Qsz,Qsinv,Qpx,Qpy,Qpz,'spline');
    Qsinv=Qsinvnew;
%     scatterQmod(Qpnx,Qpny,Qpnz,Qsinv,2);caxis([0 0.05]);
    % INTERPOLATE Qp/Qs TO Qp GRID
%     scatterQmod(Qsnx,Qsny,Qsnz,QpQs,2);caxis([0 2.5]);
    QpQsnew=interp3(Qsx,Qsy,Qsz,QpQs,Qpx,Qpy,Qpz,'spline');
    QpQs=QpQsnew;
%     scatterQmod(Qpnx,Qpny,Qpnz,QpQs,2);caxis([0 2.5]);
%     scatterQmod(Qpnx,Qpny,Qpnz,Qpinv,2);caxis([0 0.05]);
    imax=Qpnx;jmax=Qpny;kmax=Qpnz;
    Qx=Qpx;Qy=Qpy;Qz=Qpz;
else
    imax=Qsnx;jmax=Qsny;kmax=Qsnz;
    Qx=Qsx;Qy=Qsy;Qz=Qsz;
end

 
%% APPLY SMOOTHING ON V MODELS
Vpsm=smooth3(Vpmod,'gaussian',ceil(size(Vpmod,1)/size(Qpinv,1))*2+1);
Vssm=smooth3(Vsmod,'gaussian',ceil(size(Vsmod,1)/size(Qsinv,1))*2+1);
% scatterQmod(nx,ny,nz,Vpsm,2);

% INTERPOLATE V TO Q GRID
Vpnew=interp3(Vx,Vy,Vz,Vpsm,Qx,Qy,Qz,'spline');
Vsnew=interp3(Vx,Vy,Vz,Vssm,Qx,Qy,Qz,'spline');
Vratio=4/3.*Vsnew.*Vsnew./Vpnew./Vpnew;
% scatterQmod(Qsnx,Qsny,Qsnz,Vpnew,2);

%% CALCULATE Qm/Qk AND Qk
% QmQk=zeros(imax*jmax*kmax,1);
% Qkinv=zeros(imax*jmax*kmax,1);
QmQkpri=-111;Qkinvpri=-111;Qkinvstdpri=0;
QmQk=zeros(size(Qsinv));
Qkinv=zeros(size(Qsinv));
QmQkstd=zeros(size(Qsinv));
Qkinvstd=zeros(size(Qsinv));
for i=1:imax
    for j=1:jmax
        for k=1:kmax
            if (Qpinv(j,i,k)<minQ || Qsinv(j,i,k)<minQ)
                QmQk(j,i,k)=NaN;
                Qkinv(j,i,k)=NaN;
                QmQkstd(j,i,k)=NaN;
                Qkinvstd(j,i,k)=NaN;
            else
                QmQk(j,i,k)=(1-Vratio(j,i,k)*QpQs(j,i,k))/((1-Vratio(j,i,k))*QpQs(j,i,k));
                Qkinv(j,i,k)=(1-Vratio(j,i,k)*QpQs(j,i,k))/(1-Vratio(j,i,k))*Qpinv(j,i,k);
                QmQkstd(j,i,k)=NaN;
                Qkinvstd(j,i,k)=sqrt((Vratio(j,i,k)*Qpinv(j,i,k)*QpQsstd(j,i,k)/(1-Vratio(j,i,k))).^2+...
                    ((1-Vratio(j,i,k)*QpQs(j,i,k))/(1-Vratio(j,i,k))*Qpinvstd(j,i,k)).^2+...
                    ((1-QpQs(j,i,k))*Qpinv(j,i,k)*Vratio_err/((1-Vratio(j,i,k)).^2)).^2);
%                 QmQk(j,i,k)=(Qpinv(j,i,k)/Qsinv(j,i,k)-Vratio(j,i,k))/(1-Vratio(j,i,k));
%                 Qkinv(j,i,k)=(Qpinv(j,i,k)-Qsinv(j,i,k)*Vratio(j,i,k))/(1-Vratio(j,i,k));
            end
        end
    end
end
% QmQk=smooth3(QmQk);Qkinv=smooth3(Qkinv);
% scatterQmod(Qsnx,Qsny,Qsnz,QmQk,2);
% scatterQmod(imax,jmax,kmax,Qkinv,2);caxis([0 0.04]);
scatterQmod(imax,jmax,kmax,Qkinvstd,2);caxis([0 0.015]);

%% OUTPUT Qm/Qk AND Qk
tomoQmQk=zeros(imax*jmax*kmax,1);
tomoQkinv=zeros(imax*jmax*kmax,1);
tomoQkinvstd=zeros(imax*jmax*kmax,1);
for i=1:imax
    for j=1:jmax
        for k=1:kmax
            node=(i-1)*jmax*kmax+(j-1)*kmax+k;
            if isnan(QmQk(j,i,k))
                tomoQmQk(node)=QmQkpri;
            else
                tomoQmQk(node)=QmQk(j,i,k);
            end
            if isnan(Qkinv(j,i,k))
                tomoQkinv(node)=Qkinvpri;
            else
                tomoQkinv(node)=Qkinv(j,i,k);
            end
            if isnan(Qkinvstd(j,i,k))
                tomoQkinvstd(node)=Qkinvstdpri;
            else
                tomoQkinvstd(node)=Qkinvstd(j,i,k);
            end
        end
    end
end
tomoQmQk=[tomoQmQk;zeros(nQs1D,1)+1.75];
tomoQkinv=[tomoQkinv;zeros(nQs1D,1)+1.75];
tomoQkinvstd=[tomoQkinvstd;zeros(nQs1D,1)];
save(QmQkfl,'tomoQmQk','-ascii');
save(Qkfl,'tomoQkinv','-ascii');
save(Qkstdfl,'tomoQkinvstd','-ascii');