% Ploting the raypath
clear;
% workdir='/P/weisq/attentomo/3D1Dtomo/synmodv25h30';
% workdir='/Users/sowei/Work/Lau/Qtomo/3D1Dtomo/v2550h80top2/sameQpQs_fcs';
% workdir='~/GoogleDriveMSU/Work/Lau/Qtomo/3D1Dtomo/v2550h80top2_QsQps/sameQpQs_wrongerr/';
% workdir='~/GoogleDriveMSU/Work/Lau/Qtomo/3D1Dtomo/v2550h80top2_QsQps/sameQpQs/';
workdir='~/GoogleDriveMSU/Work/Lau/Qtomo/3D1Dtomo/v25h30top2_Qp/fcp0.5-20MPa0.27/';
outdir='~/GoogleDriveMSU/Work/Lau/Qtomo/fig/';
phase='p';

% Read station
staid=fopen([workdir '/cartesian.sta'], 'r');
rawsta=textscan(staid,'%5s %12.4f %12.4f %12.4f\n');
fclose(staid);
stnm=rawsta{1};
xsta=rawsta{2};
ysta=rawsta{3};
zsta=rawsta{4};

% Read raypath
rawpath=load([workdir '/raypaths.' lower(phase)]);
% rawpath=load([workdir '/raypaths.k']);
nray=0;
for i=1:size(rawpath,1)
    if (rawpath(i,1)==1234567)
        nray=nray+1;
        npt=0;
        if (upper(phase)=='S' || upper(phase)=='K')
            QpQs(nray)=rawpath(i,2);
        end
        aveQinv(nray)=rawpath(i,3);
    else
        npt=npt+1;
        x(nray,npt)=rawpath(i,1);
        y(nray,npt)=rawpath(i,2);
        z(nray,npt)=rawpath(i,3);
    end
end

% % Ploting raypath and stations
% figure(8);clf;
% set(gcf,'Renderer','Painters', 'DefaultAxesFontSize', 16,'DefaultLineLineWidth',1);
% axes('ZDir','reverse');hold on;
% if (upper(phase)=='P')
%     Qscale=1:30;
%     caxis([1 30]);
% elseif (upper(phase)=='S')
%     Qscale=1:30;
%     caxis([1 30]);
% elseif (upper(phase)=='K')
%     Qscale=1:20;
%     caxis([1 20]);
% end 
% cpt=colormap(jet(length(Qscale)));
% for i=1:nray
% % for i=3183:3183
%     if aveQinv(i)>=10
%         linew=3;
%     else
%         linew=1;
%     end
%     if (aveQinv(i)>=15 && nonzeros(z(i,1))<=700)
%         [~,idx]=min(abs(Qscale-aveQinv(i)));
%         if length(nonzeros(z(i,:)))==length(nonzeros(x(i,:)))
%         plot3(nonzeros(x(i,:)),nonzeros(y(i,:)),nonzeros(z(i,:)),...
%             'Color',cpt(idx,:),'linewidth',linew);
%         else
%             plot3(nonzeros(x(i,:)),nonzeros(y(i,:)),[nonzeros(z(i,:));0],...
%             'Color',cpt(idx,:),'linewidth',linew);
%         end
%     %     plot3(nonzeros(x1(i,:)),nonzeros(y1(i,:)),nonzeros(z1(i,:)),'m');
%         plot3(nonzeros(x(i,1)),nonzeros(y(i,1)),nonzeros(z(i,1)),'k.',...
%         'MarkerSize',30);
%     end
% end
% colorbar;
% h=colorbar;set(get(h,'title'),'string',['Path-average 1000/Q' lower(phase)]);
% plot3(xsta,ysta,zsta,'^','MarkerFaceColor','k','MarkerSize',20,'Color','k');
% xlabel('x (km)');ylabel('y (km)');zlabel('z (km)');
% view([15 5]);
% % axis([0 1000 0 400 0 700]);
% axis equal;box on;
% set(gca,'FontSize',32,'Ytick',0:400:800,'Xcolor','k','Ycolor','k','Zcolor','k');
% axis([700 1050 0 1500 0 300]);
% view([0 0]);
% % print('-depsc2',[workdir '/3Draypath' phase '.eps']);
% % saveas(gcf,[outdir '3DraypathQ' phase '.eps'],'epsc');

if (upper(phase)=='S')
    figure(8);clf;
    set(gcf,'Renderer','Painters', 'DefaultAxesFontSize', 16,'DefaultLineLineWidth',1);
%     figure('Renderer', 'Painters');clf;
    axes('ZDir','reverse');hold on;
%     box on;grid on;
    QpQsscale=0.5:.1:2.5;
    cpt=colormap(flipud(jet(length(QpQsscale))));
    caxis([0.5 2.5]);
    for i=1:nray
    % for i=3183:3183
        if QpQs(i)>=1.5
            linew=3;
        else
            linew=1;
        end
        if (QpQs(i)<=1 && aveQinv(i)>=0 && (aveQinv(i)/QpQs(i))>=0 && nonzeros(z(i,1))<=1000)
%         tmpx=nonzeros(x(i,:));
%         if (QpQs(i)<=1.5 && nonzeros(z(i,1))<=1000 && tmpx(end)<600)
            [~,idx]=min(abs(QpQsscale-QpQs(i)));
            if length(nonzeros(z(i,:)))==length(nonzeros(x(i,:)))
            plot3(nonzeros(x(i,:)),nonzeros(y(i,:)),nonzeros(z(i,:)),...
                'Color',cpt(idx,:),'linewidth',linew);
            else
                plot3(nonzeros(x(i,:)),nonzeros(y(i,:)),[nonzeros(z(i,:));0],...
                'Color',cpt(idx,:),'linewidth',linew);
            end
        %     plot3(nonzeros(x1(i,:)),nonzeros(y1(i,:)),nonzeros(z1(i,:)),'m');
            plot3(nonzeros(x(i,1)),nonzeros(y(i,1)),nonzeros(z(i,1)),'k.',...
            'MarkerSize',30);
        end
    end
    h=colorbar;set(get(h,'title'),'string','Path-average Qp/Qs');
    plot3(xsta,ysta,zsta,'^','MarkerFaceColor','k','MarkerSize',20,'Color','k');
    xlabel('x (km)');ylabel('y (km)');zlabel('z (km)');
    view([15 5]);
    axis equal;box on;
    axis([700 1050 0 1500 0 300]);view([0 0]);

    set(gca,'FontSize',36,'Ytick',0:400:800,'Xcolor','k','Ycolor','k','Zcolor','k');
% % %     saveas(gcf,'/Users/sowei/Work/Lau/Qtomo/fig/3DraypathQpQs.pdf')
% % %     print('-depsc2','/Users/sowei/Work/Lau/Qtomo/fig/3DraypathQpQs.eps');
%     saveas(gcf,[outdir '3DraypathQpQs.eps'],'epsc');
end

% 
% % Ploting the hit count
% % imax=21;jmax=9;kmax=15;dnode=50;  % nnx, nny, nnzQ horizontal 50 km, vertical 50 km
% % imax=34;jmax=14;kmax=29;  % nnx, nny, nnzQ
% % imax=21;jmax=9;kmax=21;dnode=50;  % horizontal 50 km, vertical 25-50 km
% % imax=34;jmax=14;kmax=15;dense=0;dnode=25;    % horizontal 25 km, vertical 50 km
% % imax=38;jmax=19;kmax=12;dense=2;dnode=30;zdnode=25;kmax1d=16;  % horizontal 30 km, vertical 25 km
% % imax=15;jmax=7;kmax=5;dense=0;dnode=80;zdnode=70;kmax1d=4;  % horizontal 80 km, vertical 70 km
% imax=15;jmax=7;kmax=7;dense=0;dnode=80;zdnode=50;kmax1d=4;  % horizontal 80 km, vertical 70 km
% hits=zeros(imax*jmax*kmax,1);
% xx=zeros(imax*jmax*kmax,1);
% yy=zeros(imax*jmax*kmax,1);
% zz=zeros(imax*jmax*kmax,1);
% data=load([workdir '/hits' phase]);
% for j=1:jmax
%     for i=1:imax
%         for k=1:kmax
%             node=(j-1)*imax*kmax+(i-1)*kmax+k;
%             xx(node)=i;
%             yy(node)=j;
%             zz(node)=k;
%             hits(node)=data((i-1)*jmax+j,k);
%         end
%     end
% end
% hitss=hits+.01;
% 
% figure(2);clf;
% colormap(flipud(hot));
% axes('ZDir','reverse');hold on;
% depth=zz*zdnode;
% % for i=1:length(depth)
% %     if depth(i)<=13
% %         depth(i)=(depth(i)-1)*25;
% %     else
% %         depth(i)=300+(depth(i)-13)*50;
% %     end
% % end
% % depth=(zz-1)*50;
% scatter3((xx-1)*dnode,(yy-1)*dnode,depth,hitss,hits,'filled');hold on;
% % for i=1:nray
% %     plot3(nonzeros(x(i,:)),nonzeros(y(i,:)),nonzeros(z(i,:)),'r','linewidth',3);
% % end
% axis equal;box on;
% xlabel('x (km)');ylabel('y (km)');zlabel('z (km)');
% colorbar;
% view([-30 15]);
% % axis([0 1000 0 400 0 700]);
% caxis([0 20]);
% % print('-depsc2',[workdir '/hits' phase '.eps']);
% 
% figure(4);clf;
% colormap(flipud(hot));
% depth=zz*zdnode;
% ind=find(depth==50);
% scatter((xx(ind)-1)*dnode,(yy(ind)-1)*dnode,hitss(ind),hits(ind),'filled');hold on;
% % for i=1:nray
% %     plot3(nonzeros(x(i,:)),nonzeros(y(i,:)),nonzeros(z(i,:)),'r','linewidth',3);
% % end
% axis equal;box on;
% xlabel('x (km)');ylabel('y (km)');
% colorbar;
% % axis([0 1000 0 400]);
%  caxis([0 50]);
