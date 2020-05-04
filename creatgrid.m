clear;

theta=atan(1/2);        % orientation counterclockwise from the north
% c1=[-19.0 176.75];c2=[-19.3 186.9];
% c1=[-19 176.7];c2=[-19 187.5]; % spacing 100 km
c1=[-19.0 176.5];c2=[-19.3 186.8];
% delta=0.41/cos(theta);              % grid space
% dnode=delta/180*pi*6377
dnode=60*cos(theta);                      % Grid spacing (km)
delta=dnode*180/pi/6377;

londiff=c2(2)-c1(2);
latdiff=c2(1)-c1(1);
imax=int8((londiff-latdiff*tan(theta))*cos(theta)/delta);
jmax=int8(((londiff-latdiff*tan(theta))*sin(theta)+latdiff/cos(theta))/delta);

% innlat=abs(c1(1)-c2(1));
% innlon=abs(c1(2)-c2(2));
% inna=(innlat-(innlon*tan(theta)))*cos(theta);
% innb=((innlat-(innlon*tan(theta)))*sin(theta)+innlon/cos(theta));
% jmax=int8(innb/delta+1);  % j: column number
% imax=int8(inna/delta+1);  % i: row number
% 
lat=zeros(imax,jmax);
lon=zeros(imax,jmax);

for i=1:imax
    if i==1
        lat(i,1)=c1(1);
        lon(i,1)=c1(2);
    else
        lat(i,1)=lat(i-1,1)-delta*1*sin(theta);
        lon(i,1)=lon(i-1,1)+delta*1*cos(theta);
    end
    for j=2:jmax
        lat(i,j)=lat(i,j-1)+delta*1*cos(theta);
        lon(i,j)=lon(i,j-1)+delta*1*sin(theta);
    end
end

nmax=imax*jmax;
nodes=zeros(nmax,4);
n=1;
for i=1:imax
    for j=1:jmax
        nodes(n,2)=lat(i,j);
%         if lon(i,j) > 180
%             nodes(n,3)=lon(i,j)-360;
%         else
            nodes(n,3)=lon(i,j);
%         end
        nodes(n,1)=n;
        nodes(n,4)=1;
        n=n+1;
    end
end

% Get Q nodes from Velocity nodes
j=1;
for i=1:length(nodes)
    if (mod(nodes(i,1),jmax)>=18) && (mod(nodes(i,1),jmax)<=36)
        Qnodes(j,:)=nodes(i,:);
        j=j+1;
    end
end