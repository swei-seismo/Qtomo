clear;
% Creat Qp/Qs model for checker board test
space='v2550h80';
QpQs=1;

if strcmp(space,'v70h80')
    imax=15;jmax=7;kmax=4;dnode=80;zdnode=70;  % horizontal 80 km, vertical 70 km
    width=ceil(500/zdnode);    % Width of checker board blocks (in nodes)
    depth=[60    130   200   270];
    if QpQs==1
        Qinvmax=2.25;
    else
        Qinvmax=50/1000;
    end
elseif strcmp(space,'v2550h80')
    imax=15;jmax=7;kmax=7;dnode=80;zdnode=50;  % horizontal 80 km, vertical 70 km
    width=ceil(250/zdnode);    % Width of checker board blocks (in nodes)
    depth=[25    50    75    100   150   200   300];
    if QpQs==1
        Qinvmax=2.25;
    else
        Qinvmax=50/1000;
    end
elseif strcmp(space,'v125h155')
    imax=8;jmax=4;kmax=3;dense=0;dnode=100;zdnode=100;kmax1d=4;  % horizontal 155 km, vertical 125 km
    width=ceil(250/zdnode);    % Width of checker board blocks (in nodes)
    depth=[50 175 300];
    if QpQs==1
        Qinvmax=2.25;
    else
        Qinvmax=50/1000;
    end
elseif strcmp(space,'v50100h155')
    imax=4;jmax=4;kmax=5;dnode=155;zdnode=50;  % horizontal 80 km, vertical 70 km
    width=ceil(100/zdnode);    % Width of checker board blocks (in nodes)
    depth=[50    100   150   200   300];
    if QpQs==1
        Qinvmax=2.25;
    else
        Qinvmax=50/1000;
    end
elseif strcmp(space,'v50h60')
    imax=19;jmax=9;kmax=6;dnode=60;zdnode=50;  % horizontal 80 km, vertical 70 km
    width=ceil(200/zdnode);    % Width of checker board blocks (in nodes)
    depth=[50    100   150   200   250   300];
    if QpQs==1
        Qinvmax=2.25;
    else
        Qinvmax=50/1000;
    end
elseif strcmp(space,'v2550h60s')
    imax=7;jmax=9;kmax=7;dnode=60;zdnode=50;  % horizontal 80 km, vertical 70 km
    width=ceil(150/zdnode);    % Width of checker board blocks (in nodes)
    depth=[25    50    75    100   150   200   300];
    if QpQs==1
        Qinvmax=2.25;
    else
        Qinvmax=50/1000;
    end
elseif strcmp(space,'v26h30')
    imax=38;jmax=19;kmax=11;dense=2;dnode=30;zdnode=26;  % horizontal 30 km, vertical 26 km
    width=ceil(100/zdnode);    % Width of checker board blocks (in nodes)
    depth=[38    64    90    116   142   168   194   220   246   272   298];
    Qinvmax=30/1000;
elseif strcmp(space,'v25h30')
    imax=38;jmax=19;kmax=12;dense=2;dnode=30;zdnode=25;  % horizontal 30 km, vertical 26 km
    width=ceil(100/zdnode);    % Width of checker board blocks (in nodes)
    depth=[25    50    75    100   125   150   175   200   225   250   275   300];
    Qinvmax=30/1000;
end
% nnz1D=16;   % number of layers in 1D model beneath 300 km

% % 3D part
if QpQs==1
    Qinvmin=1.;     % Qp/Qs = 1.75 or 2.25
%     outfl=['~/GoogleDriveMSU/Work/Lau/Qtomo/input/checkerQps_' space '.in'];
    outfl=sprintf('~/GoogleDriveMSU/Work/Lau/Qtomo/input/checkerQps%.2f_%s.in',Qinvmin,space);
else
    Qinvmin=0/1000;     % 1000/Q = 20 or 0
    outfl=['~/GoogleDriveMSU/Work/Lau/Qtomo/input/checker_' space '.in'];
end
% tomo=Qpri(imax*jmax*kmax,kmax,dense);       % a priori model
tomo=zeros(imax*jmax*kmax,1);
% tomo=zeros(imax*jmax*kmax,1)+Qinvmin;
for i=1:imax
    nx=floor((i)/width);
    if mod(nx,2)==0
        signx=1;
    else
        signx=-1;
    end
    for j=1:jmax
        ny=ceil((j-1)/width);
        if mod(ny,2)==0
            signy=1;
        else
            signy=-1;
        end
        for k=1:kmax
            node=(i-1)*jmax*kmax+(j-1)*kmax+k;
            nz=ceil((k)/width);
            if mod(nz,2)==0
                signz=1;
            else
                signz=-1;
            end
            if signx*signy*signz<0
                tomo(node)=Qinvmax;
            elseif signx*signy*signz>0
                tomo(node)=Qinvmin;
            end
        end
    end
end
% % Lithosphere
% for i=1:imax
%     for j=1:jmax
%         k=1;
%         node=(i-1)*jmax*kmax+(j-1)*kmax+k;
%         tomo(node)=0;
%     end
% end

% Ploting each nodes
xx=zeros(imax*jmax*kmax,1);
yy=zeros(imax*jmax*kmax,1);
zz=zeros(imax*jmax*kmax,1);
for i=1:imax
    for j=1:jmax
        for k=1:kmax
            node=(i-1)*jmax*kmax+(j-1)*kmax+k;
            xx(node)=i;
            yy(node)=j;
            zz(node)=depth(k);
        end
    end
end

figure(11);clf;
% colormap(flipud(hot));
colormap(jet);
axes('ZDir','reverse');hold on;
scatter3((xx-1)*dnode,(yy-1)*dnode,zz,100,tomo.*1000,'filled');hold on;
axis equal;box on;
xlabel('x (km)');ylabel('y (km)');zlabel('z (km)');
colorbar;
view([0 0]);
axis([0 (imax-1)*dnode 0 (jmax-1)*dnode 0 max(depth)]);
caxis([-5 Qinvmax*1000+10]);

save(outfl,'tomo','-ascii');


% % % 1D part
% outfl2='/Users/SWEI/Dropbox/Research/Lau/attentomo/input/checker1D.in';
% atten1D=zeros(kmax+nnz1D,1);
% z1D=((1:kmax+nnz1D)'-1)*50;
% for i=1:kmax+nnz1D
%     nz=floor(i/4);
%     if mod(nz,2)==0
%         signz=1;
%     else
%         signz=-1;
%     end
%     if signz>0
%         atten1D(i)=1/200;
%     elseif signz<0
%         atten1D(i)=1/1000;
%     end  
% end
% Q1D=1./atten1D;
% 
% figure(12);clf;
% plot(Q1D,z1D);
% axis ij;
% xlim([0 1200]);
% ylim([0 700]);
% 
% save(outfl2,'Q1D','-ascii');