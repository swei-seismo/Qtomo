% Invert attenuation using the maximum likelihood linear inverse method
clear;
tic;

phaselst={'P','S'};
stress={'sameQpQs_checker2'};
% stress={'fcp0.5-20MPa0.27'};
for istress=1:length(stress)
%     workdir=['/P/weisq/attentomo/3D1Dtomo/v70h80top0/' char(stress(istress))];
%     workdir=['/P/weisq/attentomo/3D1Dtomo/v2550h80top2/' char(stress(istress))];
%     workdir=['/Users/sowei/Work/Lau/Qtomo/3D1Dtomo/v125h155top2/' char(stress(istress))];
%     workdir=['/Users/sowei/Work/Lau/Qtomo/3D1Dtomo/QpQs/' char(stress(istress))];
    workdir='~/GoogleDriveMSU/Work/Lau/Qtomo/3D1Dtomo/v2550h80top2_QsQps/sameQpQs';

for iph=1:length(phaselst)
% for iph=1
    phase=char(phaselst(iph));

%% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
% % % % % % % % % % % % % Read matrices for inversion % % % % % % % % % % %
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
% imax=15;jmax=7;kmax=5;dense=0;dnode=80;zdnode=70;kmax1d=4;  % horizontal 80 km, vertical 70 km
imax=15;jmax=7;kmax=7;dense=0;dnode=80;zdnode=70;kmax1d=4;  % horizontal 80 km, vertical 25-50 km, top 2
% imax=15;jmax=7;kmax=8;dense=0;dnode=80;zdnode=70;kmax1d=4;  % horizontal 80 km, vertical 25-50 km, top 0
% imax=12;jmax=6;kmax=4;dense=0;dnode=100;zdnode=100;kmax1d=4;  % horizontal 100 km, vertical 100 km
% imax=8;jmax=4;kmax=3;dense=0;dnode=100;zdnode=100;kmax1d=4;  % horizontal 155 km, vertical 125 km
% imax=19;jmax=9;kmax=6;dense=0;dnode=60;zdnode=50;kmax1d=4;  % horizontal 80 km, vertical 70 km
ET=0.03;    % Uncertainty of theory: 0.03 s
sm1d=2;     % Smoothing factor for 1D model

% Loading matrix from Fortran codes
d=load([workdir '/Dmat' phase]);     % data
Ed=load([workdir '/ErrDmat' phase]); % variance of data
G=load([workdir '/Gmat' phase]);     % G matrix
G=sparse(G);
N=length(d);    % number of data
M=size(G,2);    % number of model parameters
% Read smoothing G matrix for 3D Q model
H3Dorig=load([workdir '/Gmat' phase 'sm']);
% % % Create smoothing G matrix for 1D Q model
H1Dorig=zeros(kmax1d,kmax1d);
H1Dorig(1,1)=1;H1Dorig(1,2)=-0.5;
H1Dorig(kmax1d,kmax1d-1)=-1;H1Dorig(kmax1d,kmax1d)=1;
for k=2:kmax1d-1
    H1Dorig(k,k-1)=-0.5;
    H1Dorig(k,k)=1;
    H1Dorig(k,k+1)=-0.5;
end
Horig=blkdiag(H3Dorig,H1Dorig);
for i=1:imax        % % Smooth between bottom layer of 3D and top layer of 1D
    for j=1:jmax
        k=kmax;
        node=(j-1)*imax*kmax+(i-1)*kmax+k;
        Horig(M-kmax1d+1,node)=-0.5/(imax*jmax);
    end
end
% Loading hits number
hitsdata=load([workdir '/hits' phase]);
hitsdata1D=load([workdir '/hits' phase '1D']);
dampratio=max(hitsdata1D)/mean(mean(hitsdata));   % ratio of damping coefficient (more damping for 1D model)
% dampratio=max(max(hitsdata))./(hitsdata1D+1);   % ratio of damping coefficient (more damping for 1D model)
hits=zeros(M,1);
for j=1:jmax
    for i=1:imax
        for k=1:kmax
            node=(j-1)*imax*kmax+(i-1)*kmax+k;
            hits(node)=hitsdata((i-1)*jmax+j,k)+0.01;
        end
    end
end
hits(imax*jmax*kmax+1:M)=hitsdata1D';

% Preparing for inversion
H1=sparse(Horig);              % smooth G
H2=sparse(blkdiag(eye(size(H3Dorig,1)),(dampratio)*eye(size(H1Dorig,1)))); % Prior G
% H2=sparse(blkdiag(eye(size(H3Dorig,1)),diag(dampratio.*dampratio))); % Prior G
% H2=sparse(eye(size(Horig,1))); % Prior G
% H2=sparse(diag(hits)); % Prior G
H=[H1;H2];
vard=Ed.^2+ET.^2;     % variance of data + uncertainty of theory (0.03s)
vard=ones(size(vard));
% vard=zeros(length(vard),1);
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
% % % % % % % % % % % % % Read matrices for inversion % % % % % % % % % % %
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %


% tic;
%% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
% % % % % % % % % % % % % Inversion with best coefficients  % % % % % % % %
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
clear ax1 ax2 ax3 ax4 hh1 hh3 hh4 para;
% attenprem=atten_prem(imax,jmax,kmax,kmax1d,2);   % a priori model = PREM
attenprem=atten_wiens08(imax,jmax,kmax,kmax1d,iph);   % a priori model = PREM
% mpri=zeros(M,1);
mpri=attenprem;
kk=0;kkk=0;

% sigmah1=10.^(-3:0.2:0);  % Smaller sigmah1 means smoother
sigmah1=8.0e-03;
% sigmah1=0.025:0.005:0.035;
% sigmah1=[1e-3 0.005:0.005:0.02 0.04:0.02:0.1];
% sigmah1=2.5e-2:5e-3:4.5e-2;
expo2=[4];     % Larger expo2 and sigmah2 means less damping
reserr=zeros(length(sigmah1)*length(expo2),4);
ireserr=1;
for ii1=1:length(sigmah1)   % Different smoothing
    varh1=(sigmah1(ii1)^2)*ones(size(H1,2),1);
%     mpri=Qpri(imax*jmax*kmax,kmax,dense);       % a priori model
    for ii2=1:length(expo2)
%         mpri=zeros(M,1);       % a priori model = no attenuation
        mpri=attenprem;
        outfl=sprintf('%s/Qinv_%.1e_%.0e.%s',workdir,sigmah1(ii1),10^(expo2(ii2)),lower(phase));
%         if exist(outfl,'file')
%             continue
%         end
        sprintf('%.1e_%.0e',sigmah1(ii1),10^(expo2(ii2)))
        sigmah2=10^expo2(ii2);   % Different damping
        varh2=(sigmah2^2)*ones(size(H2,2),1);
%         for inode=1:(imax*jmax)     % Fix the lithosphere
%             ind1=1+(inode-1)*kmax;
%             varh2(ind1)=varh2(ind1)*1e-4;
% %             mpri(ind1)=attenprem(1);
%         end
        varh=[varh1;varh2];

        clear GG dobs;
%         G=sparse(G);
%         H=sparse(H);
        Cd=sparse(diag(vard));
        Ch=sparse(diag(varh));
        if condest(Ch)>1e15
            continue
        end
        GG = G'*(Cd\G) + H'*(Ch\H);
        GG = sparse(GG);
        dobs = G'*(Cd\d) + H'*(Ch\H)*mpri;
%         options = optimset('MaxIter',1000*length(dobs),'TolX',100);
        options = optimset('MaxIter',1000*N);
        [mest,resnorm,residual,exitflag]=lsqnonneg(GG,dobs,options);
        if exitflag==0
            sprintf('skip');
            continue
        end

% % %   Model errors
        Ginv=inv(GG'*GG);
        merr=sqrt(resnorm/(length(d)-length(mest)-1)*diag(Ginv));
        merr=full(merr);
        
        % Output
%         save(outfl,'mest','-ascii');
        save([outfl '_varm'],'merr','-ascii');
    end
end

% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
% % % % % % % % % % % % % Inversion with best coefficients  % % % % % % % %
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
end
end
tend=toc;
fprintf('%d minutes and %f seconds\n',floor(tend/60),rem(tend,60));