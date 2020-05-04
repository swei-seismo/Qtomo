% % % Jointly invert for Qp/Qs based on 1/Qp and t*(S)
clear;
tic;
% workdir='/P/weisq/attentomo/3D1Dtomo/v150h155/fcs0.5-20MPa0.27_svd';

% stress={'Qps_checkerclean'};
stress={'sameQps_checkersvd_vard2.25'};
% stress={'Qps'};
for istress=1:length(stress)
%% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
% % % % % % % % % % % % % Read matrices for inversion % % % % % % % % % % %
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
%     workdir=['/P/weisq/attentomo/3D1Dtomo/v70h80top0/' char(stress(istress))];
%     workdir=['/Users/sowei/Work/Lau/Qtomo/3D1Dtomo/QpQs/same' char(stress(istress))];
    workdir=['~/GoogleDriveMSU/Work/Lau/Qtomo/3D1Dtomo/v2550h80top2_QsQps/' char(stress(istress))];

% imax=8;jmax=4;kmax=3;dense=0;dnode=100;zdnode=100;kmax1d=4;  % horizontal 155 km, vertical 125 km
imax=15;jmax=7;kmax=7;dense=0;dnode=80;zdnode=70;kmax1d=4;  % horizontal 80 km, vertical 25-50 km
% imax=15;jmax=7;kmax=8;dense=0;dnode=80;zdnode=70;kmax1d=4;  % horizontal 80 km, vertical 25-50 km, top 0
% imax=21;jmax=9;kmax=15;dense=0;  % nnx, nny, nnzQ
% imax=34;jmax=14;kmax=29;dense=1;  % horizontal 25 km, vertical 25 km
% imax=19;jmax=9;kmax=6;dense=0;dnode=60;zdnode=50;kmax1d=4;  % horizontal 60 km, vertical 50 km
% imax=15;jmax=7;kmax=4;dense=0;dnode=80;zdnode=70;kmax1d=4;  % horizontal 80 km, vertical 70 km
% imax=8;jmax=4;kmax=3;dense=0;dnode=50;zdnode=150;kmax1d=3;  % horizontal 155 km, vertical 150 km
% imax=15;jmax=7;kmax=5;dense=0;dnode=80;zdnode=70;kmax1d=4;  % horizontal 80 km, vertical 70 km
ET=0.03;    % Uncertainty of theory: 0.03 s

Qpinv=load([workdir '/Qinv.p']);
d=load([workdir '/DmatPS']);
Ed=load([workdir '/ErrDmatPS']); % variance of data
vard=Ed.^2+ET.^2;     % variance of data + uncertainty of theory (0.03s)
G=load([workdir '/GmatPS']);
N=length(d);    % number of data
M=size(G,2);    % number of model parameters
WG=G;
for ig=1:N
    WG(ig,:)=WG(ig,:)./vard(ig);
end
% G=sparse(G);
% WG=sparse(WG);
% wd=load([workdir '/WDmatSP']);
smd=load([workdir '/DmatPSsm']);
% smd=[smd;zeros(kmax1d,1)];
% smd=sparse(smd);
smG=load([workdir '/GmatPSsm']);
smG1D=zeros(kmax1d,kmax1d);
smG1D(1,1)=1;smG1D(1,2)=-0.5;
smG1D(kmax1d,kmax1d-1)=-1;smG1D(kmax1d,kmax1d)=1;
for k=2:kmax1d-1
    smG1D(k,k-1)=-0.5;
    smG1D(k,k)=1;
    smG1D(k,k+1)=-0.5;
end
smG=blkdiag(smG,smG1D);
for i=1:imax        % % Smooth between bottom layer of 3D and top layer of 1D
    for j=1:jmax
        k=kmax;
        node=(j-1)*imax*kmax+(i-1)*kmax+k;
        smG(M-kmax1d+1,node)=-0.5/(imax*jmax);
    end
end
% smG=sparse(smG);
% smwd=load([workdir '/WDmatSPsm']);
hitsdata=load([workdir '/hitsS']);
hitsdata1D=load([workdir '/hitsS1D']);
dampratio=max(hitsdata1D)/mean(mean(hitsdata));   % ratio of damping coefficient (more damping for 1D model)
hits=zeros(imax*jmax*kmax,1);
for j=1:jmax
    for i=1:imax
        for k=1:kmax
            node=(i-1)*jmax*kmax+(j-1)*kmax+k;
            hits(node)=hitsdata((i-1)*jmax+j,k);
        end
    end
end
% I=diag(ones(size(G,2),1));
I=eye(size(G,2));

% m0=dampQ(G,kmax,dense);       % a priori model
attenprem=zeros(imax*jmax*kmax+kmax1d,1)+2.25;   % a priori model

% Regulartions
% lamda: Smoothing coefficient squared, larger lamda means smooth
% epsilon: Damping coefficient squared, larger epsilon means larger damping
lamda=[10];
% lamda=0.1;
epsilon=[2 100:100:200];  
pvalues=0.01:0.02:0.09;
% % % tic;
reserr=zeros(length(lamda)*length(epsilon),6);
ireserr=1;
clear ax1 ax2 ax3 ax4 hh1 hh3 hh4 para;
ifig=1;
mpri=attenprem;
for ii1=1:length(lamda)
%     lamda=1/sigmah1(ii1);
%     for ii2=1:length(epsilon)
    for ii3=1:length(pvalues)
        mpri=attenprem;
%         basenm=sprintf('%s/Qinv_%.1f_%.1f',workdir,lamda(ii1),epsilon(ii2));
        basenm=sprintf('%s/Qinv_%.1f_p%.2f',workdir,lamda(ii1),pvalues(ii3));
        clear GG dd;
        GG=[WG;lamda(ii1)*smG];
%         GG=sparse(GG);
        dd=[d./vard-WG*mpri;zeros(size(smG,1),1)];
%         dd=sparse(dd);
        [U,S,V]=svd(GG,0);
        
% % % % %  % Damped SVD (Menke, 2014, Eq. 7.45)
%         p=size(GG,2);
%         Sp=S(1:p,1:p);
%         Vp=V(:,1:p);
%         Up=U(:,1:p);
%         damping=eye(size(S,1))*epsilon(ii2);
% %         for inode=1:(imax*jmax)     % Fix the lithosphere
% %             ind1=1+(inode-1)*kmax;
% %             epsilon(ind1,ind1)=epsilon(ind1,ind1)*1e+15;
% % %             mpri(ind1)=attenprem(1);
% %         end
%         for inode=(imax*jmax*kmax+1):M      % Fix the 1D part
%             damping(inode,inode)=damping(inode,inode)*1e+10;
%         end
%         Sp=(S(1:p,1:p))^2+damping(1:p,1:p);  % Damped SVD (Menke, 2014, Eq. 7.45)
% %         dmest=Vp/Sp*Up'*dd;
%         Gg=Vp*(Sp\(S(1:p,1:p)))*Up';
% %         dmest=Vp*(Sp\(S(1:p,1:p)))*Up'*dd;
% % % % %  % Damped SVD (Menke, 2014, Eq. 7.45)

% % % %  % SVD (Menke, 2014, Eq. 7.40)
        p=floor(size(GG,2)*pvalues(ii3));
%         p=floor(size(GG,2)/2);
%         p=size(GG,2);
        Vp=V(:,1:p);
        Up=U(:,1:p);
        Sp=S(1:p,1:p);
        Gg=Vp*(Sp\Up');
% % % %  % SVD (Menke, 2014, Eq. 7.40)

        figure(2);clf;grid on;box on;hold on;
        plot(diag(S),'bo');
        plot(p,S(p,p),'r*','markersize',50,'linewidth',5);
        xlabel('Index');ylabel('Singular value');
        set(gca,'fontsize',24);

        dmest=Gg*dd;
%         dmest=Vp*(Sp\(S(1:p,1:p)))*Up'*dd;
        mest=mpri+dmest;
%         mest=dmest;
        R=Vp*Vp';  % Model resolution
        spreadR=(norm(R-eye(size(R,1)))).^2;
        covd=diag([vard;zeros(size(smG,1),1)]);
%         Gg=Vp/Sp*Up';
        covm=Gg*covd*Gg';     % Model covariance
%         covm=Gg*Gg';     % Model covariance
        varm=diag(covm);varm=sqrt(varm);

        % Output
        outfl=[basenm '.ps'];
%         outfl
        save(outfl,'mest','-ascii');
        outflres=[outfl '_varm']; % Model variance
        save(outflres,'varm','-ascii');
        % Resolution length & Predication error 
        Gm= G*mest;
        Err=(d-Gm)'*(d-Gm);
        L=dmest'*dmest;  % Resolution length
%         reserr(ireserr,:)=[lamda(ii1),epsilon(ii2),Err,L,mean(varm),spreadR];
        reserr(ireserr,:)=[lamda(ii1),pvalues(ii3),Err,L,mean(varm),spreadR];
        ireserr=ireserr+1;
        
        % Output 1/Qs
        Qsinv=Qpinv.*mest;
        outfls=[basenm '.s'];
        save(outfls,'Qsinv','-ascii');
    end
end
% Store Q inversion misfit for each set of t*
misfl=[workdir '/Qmisfit.all'];
fid=fopen(misfl,'w');
for ii=1:size(reserr,1)
    fprintf(fid,'%10.3f %10.3f %10.3f %10.5f %10.3f %10.3e\n',reserr(ii,:));
end
fclose(fid);
end
tend=toc;
fprintf('%d minutes and %f seconds\n',floor(tend/60),rem(tend,60));

% % Output
% save([workdir '/attenP.lsqr'],'m','-ascii');
