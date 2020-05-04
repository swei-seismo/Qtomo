% Invert attenuation using the maximum likelihood linear inverse method
clear all;
global G H vard varh;
tic;

% stress={'0.5-20MPa0','0.5-20MPa027','0.5-20MPa06',...
%     '10MPa0.0','10MPa0.27','10MPa0.6',...
%     '20MPa0.0','20MPa0.27','20MPa0.6',...
%     '5MPa0.0','5MPa0.27','5MPa0.6'};
% stress={'fcp0.5-20MPa0.27','checker','synmod'};
stress={'fcp0.5-20MPa0.27'};
for istress=1:length(stress)
%% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
% % % % % % % % % % % % % Read matrices for inversion % % % % % % % % % % %
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
%     workdir=['/P/weisq/attentomo/3D1Dtomo/v25h30/' char(stress(istress))];
    workdir=['/Users/sowei/GoogleDriveMSU/Work/Lau/Qtomo/3D1Dtomo/v25h30top2_Qp/' char(stress(istress))];

% workdir='/P/weisq/attentomo/3D1Dtomo/v25h30';


% imax=34;jmax=14;kmax=29;dense=1;dnode=25;  % horizontal 25 km, vertical 25 km
% imax=21;jmax=9;kmax=21;dense=2;dnode=50;  % horizontal 50 km, vertical 25-50 km
% imax=34;jmax=14;kmax=21;dense=2;dnode=25;  % horizontal 25 km, vertical 25-50 km
% imax=34;jmax=14;kmax=15;dense=4;dnode=25;    % horizontal 25 km, vertical 50 km
% imax=21;jmax=9;kmax=15;dense=0;dnode=50;    % horizontal 50 km, vertical 50 km
imax=38;jmax=19;kmax=12;dense=2;dnode=30;zdnode=25;kmax1d=8;  % horizontal 30 km, vertical 25 km
ET=0.03;    % Uncertainty of theory: 0.03 s
sm1d=2;     % Smoothing factor for 1D model

% Loading matrix from Fortran codes
d=load([workdir '/DmatP']);     % data
Ed=load([workdir '/ErrDmatP']); % variance of data
G=load([workdir '/GmatP']);     % G matrix
G=sparse(G);
N=length(d);    % number of data
M=size(G,2);    % number of model parameters
% Read smoothing G matrix for 3D Q model
H3Dorig=load([workdir '/GmatPsm']);
% % Create smoothing G matrix for 1D Q model
H1Dorig=zeros(kmax1d,kmax1d);
H1Dorig(1,1)=sm1d;H1Dorig(1,2)=-sm1d/2;
H1Dorig(kmax1d,kmax1d-1)=-sm1d;H1Dorig(kmax1d,kmax1d)=sm1d;
for k=2:kmax1d-1
    H1Dorig(k,k-1)=-sm1d/2;
    H1Dorig(k,k)=sm1d;
    H1Dorig(k,k+1)=-sm1d/2;
end
Horig=blkdiag(H3Dorig,H1Dorig);
for i=1:imax        % % Smooth between bottom layer of 3D and top layer of 1D
    for j=1:jmax
        k=kmax;
        node=(j-1)*imax*kmax+(i-1)*kmax+k;
        Horig(M-kmax1d+1,node)=-sm1d/2/(imax*jmax);
    end
end
% Loading hits number
hitsdata=load([workdir '/hitsP']);
hitsdata1D=load([workdir '/hitsP1D']);
dampratio=max(hitsdata1D)/mean(mean(hitsdata));   % ratio of damping coefficient (more damping for 1D model)
% dampratio=max(max(hitsdata))./(hitsdata1D+1);   % ratio of damping coefficient (more damping for 1D model)
hits=zeros(M,1);
for j=1:jmax
    for i=1:imax
        for k=1:kmax
            node=(j-1)*imax*kmax+(i-1)*kmax+k;
            hits(node)=hitsdata((i-1)*jmax+j,k);
        end
    end
end
hits(imax*jmax*kmax+1:M)=hitsdata1D';
% % hits=hits.^2;

% Preparing for inversion
H1=sparse(Horig);              % smooth G
H2=sparse(blkdiag(eye(size(H3Dorig,1)),(dampratio^2)*eye(size(H1Dorig,1)))); % Prior G
% H2=sparse(blkdiag(eye(size(H3Dorig,1)),diag(dampratio.*dampratio))); % Prior G
% H2=sparse(diag(hits)); % Prior G
H=[H1;H2];
vard=Ed.^2+ET.^2;     % variance of data + uncertainty of theory (0.03s)
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
% % % % % % % % % % % % % Read matrices for inversion % % % % % % % % % % %
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %

% % tstart=tic;
%% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
% % % % % % % % % % % % % Inversion with various coefficients % % % % % % %
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
% sigmah1=10.^(-8:-1);
% % sigmah1=[0.1 0.01 0.001 0.0005 0.0001 0.00001 0.000001];  % Smaller sigmah1 means smoother
% expo2=2:8;     % Larger expo2 and sigmah2 means less damping
% reserr=zeros(length(sigmah1)*length(expo2),3);
% ireserr=1;
% for ii1=1:length(sigmah1)   % Different smoothing
%     varh1=(sigmah1(ii1)^2)*ones(size(H1,2),1);
% %     mpri=Qpri(M,kmax,dense);       % a priori model
%     for ii2=1:length(expo2)
%         mpri=zeros(M,1);       % a priori model = no attenuation
% %         outfl=[workdir '/Qinv_' num2str(sigmah1(ii1)) '_' num2str(10^(expo2(ii2))) '.p'];
% %         for sigmah2=10.^([-5:expo2 expo2-1:-1:-5])
% %         for sigmah2=10.^(expo2(ii2):-1:-5)
%         for sigmah2=10^expo2(ii2)   % Different damping
%             varh2=(sigmah2^2)*ones(size(H2,2),1);
%             for inode=1:(imax*jmax)     % Fix the lithosphere
%                 ind1=1+(inode-1)*kmax;
%                 varh2(ind1)=1e-10*(sigmah2^2);
%             end
% %             for i=1:(imax*jmax*kmae-10e-10x)    % Fix the slab
% %                 if mpri(i)==0
% %                     varh2(i)=0.08*(sigmah2^2);
% %                 end
% %             end
%             varh=[varh1;varh2];
% 
%             dobs=d-G*mpri;
%             % dmest=bicg( @glsfcn, G'*(dobs./vard), 1e-5, 4*length(dobs));
%             dmest=lsqr( @glsfcn, G'*(dobs./vard), 1e-5, 10*length(dobs));
%             mest=mpri+dmest;
%             mpri=(mest+abs(mest)+0.000001)/2;     % damp negative ones to 0.0000005
% 
%             % Ploting each nodes
%             tomo=zeros(imax*jmax*kmax,1);
%             xx=zeros(imax*jmax*kmax,1);
%             yy=zeros(imax*jmax*kmax,1);
%             zz=zeros(imax*jmax*kmax,1);
%             for i=1:imax
%                 for j=1:jmax
%                     for k=1:kmax
%                         node=(i-1)*jmax*kmax+(j-1)*kmax+k;
%                         xx(node)=i;
%                         yy(node)=j;
%                         zz(node)=k;
%                         tomo(node)=mest((i-1)*jmax*kmax+(j-1)*kmax+k);
%                     end
%                 end
%             end
% 
% %             figure(1);clf;
% %             colormap(flipud(hot));
% %         %     colormap(jet);
% %             axes('ZDir','reverse');hold on;
% %             scatter3((xx-1)*dnode,(yy-1)*dnode,(zz-1)*zdnode,20,tomo.*1000,'filled');hold on;
% %             axis equal;box on;
% %             xlabel('x (km)');ylabel('y (km)');zlabel('z (km)');
% %             colorbar;
% %             view([-30 15]);axis([0 (imax-1)*dnode 0 (jmax-1)*dnode 0 (kmax-1)*zdnode]);
% %         %     caxis([0 30]);
% 
%         end
% %         print('-depsc2',[workdir '/Qpinv_' num2str(sigmah1(ii1)) '_' num2str(10^(expo2(ii2))) '.eps']);
% %         % Output
% %         save(outfl,'mest','-ascii');
% 
% %         % size of model covariance and spread of model resolution
% %         covmest=0;
% %         spdmres=0;
% %         for k = 1:floor(M/1000):M
% %             s = zeros(M,1);
% %             s(k)=1;
% %             % compute k-th column of inverse of A=(G'*Cdi*G + H'*Chi*H)
% %             % note inv(A) symmetric, so column is also row
% %             ai = lsqr( @glsfcn, s, 1e-5, 4*N );
% %             % k-th row of the resolution matrix
% %             r_row=((ai')*(G'))./(vard')*G;
% %             % size of model covariance: sum of k-th element of k-th column
% %             % of inv(A)
% %             covmest = covmest+ai(k);
% %             % spread of model resolution
% %             spreadR_row = 0;
% %             for j = 1:M
% %                 if j == k
% %                     Rij = r_row(j)-1;
% %                 else
% %                     Rij = r_row(j);
% %                 end
% %                 spreadR_row = spreadR_row+Rij^2;
% %             end
% %             spdmres = spdmres+spreadR_row;
% %         end
%         
%         % Resolution length & Predication error 
%         Gm=G*mest;
%         Err=(d-Gm)'*(d-Gm);
% %         reserr(ireserr,:)=[sigmah1(ii1),10^expo2(ii2),Err,spdmres,covmest];
%         reserr(ireserr,:)=[sigmah1(ii1),10^expo2(ii2),Err];
%         ireserr=ireserr+1;
%     end
% end
% % tend=toc;
% % fprintf('%d minutes and %f seconds\n',floor(tend/60),rem(tend,60));
% 
% figure(4);clf;
% h=plot3(reserr(:,1),reserr(:,2),reserr(:,3),'*');hold on;
% grid on;box on;
% set(gca,'XScale','log','YScale','log');
% xlabel('\sigma_1 (1/smoothing coeff)');ylabel('\sigma_2 (1/damping coeff)');
% zlabel('Prediction error');
% title('Qp inversion with different parameters');
% saveas(h,[char(stress(istress)) '.fig']);
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
% % % % % % % % % % % % % Inversion with various coefficients % % % % % % %
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %


%% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 
% % % % % % % % % % % % % Inversion with best coefficients  % % % % % % % %
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
clear ax1 ax2 ax3 ax4 hh1 hh3 hh4 para;
% attenprem=atten_prem(imax,jmax,kmax,kmax1d,1);   % a priori model = PREM
attenprem=atten_wiens08(imax,jmax,kmax,kmax1d,1);   % a priori model = PREM
% mpri=zeros(M,1);
mpri=attenprem;
kk=0;kkk=0;
keynode=530;
keynode2=580;
% figure(7);clf;ax1=gca;
% figure(8);clf;ax2=gca;
% figure(9);clf;ax3=gca;
% figure(10);clf;ax4=gca;
% plot(ax3,mpri(end-7:end)*1000,(1:8)*50+300,'r','linewidth',10);hold(ax3,'on');
% plot(ax1,mpri((keynode-1)*kmax+1:keynode*kmax)*1000,...
%     (1:kmax)*25,'r','linewidth',10);hold(ax1,'on');
% plot(ax4,mpri((keynode2-1)*kmax+1:keynode2*kmax)*1000,...
%     (1:kmax)*25,'r','linewidth',10);hold(ax4,'on');

% sigmah1=10.^(-5:1:-1);  % Smaller sigmah1 means smoother
% sigmah1=[1e-6 1e-5 1e-4:1e-4:1e-3 1e-2 1e-1];
sigmah1=[7*1e-4];
% expo2=-6:2:-4;     % Larger expo2 and sigmah2 means less damping, no change if expo2>0
% sigmah1=(2:9).*1e-4;  % Smaller sigmah1 means smoother
expo2=[6];     % Larger expo2 and sigmah2 means less damping, no change if expo2>0
reserr=zeros(length(sigmah1)*length(expo2),4);
ireserr=1;
for ii1=1:length(sigmah1)   % Different smoothing
    varh1=(sigmah1(ii1)^2)*ones(size(H1,2),1);
%     mpri=Qpri(imax*jmax*kmax,kmax,dense);       % a priori model
    for ii2=1:length(expo2)
%         ii2=ii1;
%         tic;
%         mpri=zeros(M,1);       % a priori model = no attenuation
        mpri=attenprem;
%         outfl=[workdir '/Qinv_' num2str(sigmah1(ii1)) '_' num2str(10^(expo2(ii2))) '.p'];
        outfl=sprintf('%s/Qinv_%.1e_%.0ehigh.p',workdir,sigmah1(ii1),10^(expo2(ii2)));
        outfl1=sprintf('%s/Qinv_%.1e_%.0e.p1',workdir,sigmah1(ii1),10^(expo2(ii2)));
%         for sigmah2=10.^([-5:expo2(ii2) expo2(ii2)-1:-1:-5])
        dampiter=0;
        for sigmah2=10.^(expo2(ii2):-1:-10)
%         for sigmah2=10^expo2(ii2)   % Different damping
            dampiter=dampiter+1;
            varh2=(sigmah2^2)*ones(size(H2,2),1);
%             for inode=1:(imax*jmax)     % Fix the lithosphere
%                 ind1=1+(inode-1)*kmax;
%                 varh2(ind1)=varh2(ind1)*1e-4;
% %                 mpri(ind1)=attenprem(1);
%             end
            varh=[varh1;varh2];

            dobs=d-G*mpri;
            G=sparse(G);
            dobs=sparse(dobs);
            H=sparse(H);
            vard=sparse(vard);
            varh=sparse(varh);
            dmest=bicg( @glsfcn, G'*(dobs./vard), 1e-5, 4*length(dobs));
%             dmest=lsqr( @glsfcn, G'*(dobs./vard), 1e-5, 4*length(dobs));
%             mest=mpri+dmest;
%             dmest=(dmest+abs(dmest)+0.001)/2;     % damp negative ones to 0.0005
%             mpri=mpri+dmest;
            mest=mpri+dmest;
%             mest=dmest;
            mpri=(mest+abs(mest)+0.000001)/2;     % damp negative ones to 0.0005
%             mpri(end-kmax1d+1:end)=attenprem(end-kmax1d+1:end);
% %             if dampiter==1
% %                 save(outfl1,'mest','-ascii');
% %             end

            kk=kk+1;
%             % Ploting 1D model
%             para(kk)={sprintf('%d: %.1e %.1e',kk,sigmah1(ii1),sigmah2)};
%             if kk==1
%                 hh3(kk)=plot(ax3,mest(end-7:end)*1000,(1:8)*50+300,'linewidth',5);hold(ax3,'all');
% %                 text(ax3,mest(end-4)*1000,500,num2str(kk));hold(ax3,'all');
%                 hh1(kk)=plot(ax1,mest((keynode-1)*kmax+1:keynode*kmax)*1000,...
%                     (1:kmax)*25,'linewidth',5);hold(ax1,'all');
%                 hh4(kk)=plot(ax4,mest((keynode2-1)*kmax+1:keynode2*kmax)*1000,...
%                     (1:kmax)*25,'linewidth',5);hold(ax4,'all');
% %                 text(ax1,mest((keynode-1)*kmax+3)*1000,200,num2str(kk));hold(ax3,'all');
%             else
%                 hh3(kk)=plot(ax3,mest(end-7:end)*1000,(1:8)*50+300,'linewidth',2);hold(ax3,'all');
% %                 text(ax3,mest(end-4)*1000,500,num2str(kk));hold(ax3,'all');
%                 hh1(kk)=plot(ax1,mest((keynode-1)*kmax+1:keynode*kmax)*1000,...
%                     (1:kmax)*25,'linewidth',2);hold(ax1,'all');
%                 hh4(kk)=plot(ax4,mest((keynode2-1)*kmax+1:keynode2*kmax)*1000,...
%                     (1:kmax)*25,'linewidth',2);hold(ax4,'all');
% %                 text(ax1,mest((keynode-1)*kmax+3)*1000,200,num2str(kk));hold(ax3,'all');
%             end
% 
%             % Ploting misfit
%             Gm=G*(mpri+dmest);
%             Err=(d-Gm)'*(d-Gm);
%             plot(ax2,kk,Err,'b*');hold(ax2,'on');
            
%             % Ploting each nodes
%             tomo=zeros(imax*jmax*kmax,1);
%             xx=zeros(imax*jmax*kmax,1);
%             yy=zeros(imax*jmax*kmax,1);
%             zz=zeros(imax*jmax*kmax,1);
%             for i=1:imax
%                 for j=1:jmax
%                     for k=1:kmax
%                         node=(i-1)*jmax*kmax+(j-1)*kmax+k;
%                         xx(node)=i;
%                         yy(node)=j;
%                         zz(node)=k;
%                         tomo(node)=mest((i-1)*jmax*kmax+(j-1)*kmax+k);
%                     end
%                 end
%             end
% 
%             figure(1);clf;
% %             set(gcf,'Visible','off');
%             colormap(flipud(hot));
%         %     colormap(jet);
%             axes('ZDir','reverse');hold on;
%             scatter3((xx-1)*dnode,(yy-1)*dnode,(zz-1)*zdnode,20,tomo.*1000,'filled');hold on;
%             axis equal;box on;
%             xlabel('x (km)');ylabel('y (km)');zlabel('z (km)');
%             colorbar;
%             view([-30 15]);axis([0 (imax-1)*dnode 0 (jmax-1)*dnode 0 (kmax-1)*zdnode]);
%             caxis([0 30]);

        end
% %         print('-depsc2',[workdir '/Qpinv_' num2str(sigmah1(ii1)) '_' num2str(10^(expo2(ii2))) '.eps']);
%         figfl=sprintf('%s/Qpinv_%.0e_%.0e.eps',workdir,sigmah1(ii1),10^(expo2(ii2)));
%         print('-depsc2',figfl);
        kkk=kkk+1;
%         % Ploting 1D model
%         para(kkk)={sprintf('%d: %.1e %.1e',kkk,sigmah1(ii1),10^expo2(ii2))};
%         if kkk==1
%             hh3(kkk)=plot(ax3,mest(end-7:end)*1000,(1:8)*50+300,'linewidth',5);hold(ax3,'all');
%             hh1(kkk)=plot(ax1,mest((keynode-1)*kmax+1:keynode*kmax)*1000,...
%                 (1:kmax)*25,'linewidth',5);hold(ax1,'all');
%             hh4(kkk)=plot(ax4,mest((keynode2-1)*kmax+1:keynode2*kmax)*1000,...
%                 (1:kmax)*25,'linewidth',5);hold(ax4,'all');
%         else
%             hh3(kkk)=plot(ax3,mest(end-7:end)*1000,(1:8)*50+300,'linewidth',2);hold(ax3,'all');
%             hh1(kkk)=plot(ax1,mest((keynode-1)*kmax+1:keynode*kmax)*1000,...
%                 (1:kmax)*25,'linewidth',2);hold(ax1,'all');
%             hh4(kkk)=plot(ax4,mest((keynode2-1)*kmax+1:keynode2*kmax)*1000,...
%                 (1:kmax)*25,'linewidth',2);hold(ax4,'all');
%         end
        % Output
        save(outfl,'mest','-ascii');
%         % Resolution length & Predication error 
%         Len=mest'*mest;
        Gm=G*mest;
        Err=(d-Gm)'*(d-Gm);
        reserr(ireserr,:)=[sigmah1(ii1),10^expo2(ii2),Err,Err/N];
        ireserr=ireserr+1;
    end
end
% legend(ax1,hh1(:),para(:),'Location','East');
% % plot(ax1,mest((keynode-1)*kmax+1:keynode*kmax)*1000,(1:kmax)*25,'k','linewidth',2);
% axis(ax1,'ij');
% ylabel(ax1,'Depth (km)');xlabel(ax1,'1000/Q');
% title(ax1,['Node ' num2str(keynode) ' ' num2str(hits((keynode-1)*kmax+1))]);
% ylim(ax1,[0 300]);
% %xlim(ax1,[-5 40]);
% xlabel(ax2,'Iteration');ylabel(ax2,'Solution length');
% legend(ax3,hh3(:),para(:),'Location','East');
% % plot(ax3,mest(end-7:end)*1000,(1:8)*50+300,'k','linewidth',5);
% axis(ax3,'ij');
% ylabel(ax3,'Depth (km)');xlabel(ax1,'1000/Q');
% title(ax3,'1D beneath 300 km');
% ylim(ax3,[300 700]);
% % xlim(ax3,[-5 20]);
% legend(ax4,hh1(:),para(:),'Location','East');
% % plot(ax4,mest((keynode2-1)*kmax+1:keynode2*kmax)*1000,(1:kmax)*25,'k','linewidth',2);
% axis(ax4,'ij');
% ylabel(ax4,'Depth (km)');xlabel(ax1,'1000/Q');
% title(ax4,['Node ' num2str(keynode2) ' ' num2str(hits((keynode2-1)*kmax+1))]);
% ylim(ax4,[0 300]);
% %xlim(ax4,[-5 40]);

atten1D=mest(imax*jmax*kmax+1:M);

% Store Q inversion misfit for each set of t*
misfl=[workdir '/Qmisfit.all'];
fid=fopen(misfl,'a');
for ii=1:size(reserr,1)
    fprintf(fid,'%e %e %e %e\n',reserr(ii,:));
end
fclose(fid);
% save(misfl,'reserr','-ascii');
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
% % % % % % % % % % % % % Inversion with best coefficients  % % % % % % % %
%%% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %

end
tend=toc;
fprintf('%d minutes and %f seconds\n',floor(tend/60),rem(tend,60));
