% FDTD code in 2D case. 
% Free surface is applied at all sides
% PML is applied at left side but it is disabled by commenting
% Code last modified on 28th December 2008, 8:00 PM
% Grid structure and refence paper is given in the same folder

%% Clear work space
clc;
clear *;

%% Define simulation parameters
rho=2667;       %Properties of pure aluminium ref:http://www.sciencedirect.com/science?_ob=ArticleURL&_udi=B6TW0-3X3K03P-B&_user=518931&_rdoc=1&_fmt=&_orig=search&_sort=d&view=c&_acct=C000025838&_version=1&_urlVersion=0&_userid=518931&md5=358cd234cabafb4ee9a4f8964694194d
cl=6396;
cs=3103;
mu=cs^2*rho;
lambda=cl^2*rho-2*mu;
freq=2.25*1e6;  %Hz
wl=cl/freq; %wavelength
h=wl/30;
dt=h/cl/(1.5);

depth=20;
width=20;
M=round(depth*1e-3/h);
N=round(width*1e-3/h);
tmax=round(10e-6/dt);

%% Generate input pulse ...
% The function gethanningpulse returns 3 cycles of hanning windowed
% gaussian pulse at a given frequency
ip=gethanningpulse(freq,dt,tmax);

%% Create figure to image plot the resultant velocity matrix
uaxis=axes('parent',figure('Name','U'));
clrmap=load('rainbowcolormap.txt');

%% Few constants 
dtbyrhodh=dt/(rho*h);
lambdaplus2mudtbydh=(lambda+2*mu)*dt/h;
lambdadtbydh=lambda*dt/h;
dtmubydh=dt*mu/h;
freesurfaceconstant=lambda/(lambda+2*mu);

%% Set all the field variable to zero in the begining
u=zeros(M+2,N+2);
v=zeros(M+2,N+2);
txx=zeros(M+2,N+2);
tzz=zeros(M+2,N+2);
txz=zeros(M+1,N+1);
for tind=1:tmax
    %% FDTD Operation
    u(2:M+1,1:N+1)=u(2:M+1,1:N+1)+dtbyrhodh*(txx(2:M+1,2:N+2)-txx(2:M+1,1:N+1)+txz(2:M+1,1:N+1)-txz(1:M,1:N+1));
    v(1:M+1,2:N+1)=v(1:M+1,2:N+1)+dtbyrhodh*(txz(1:M+1,2:N+1)-txz(1:M+1,1:N)+tzz(2:M+2,2:N+1)-tzz(1:M+1,2:N+1));

   
    %% Applying input pulse
    if ip(tind)~=123456
        v(round(M/2),round(N/2))=ip(tind);
    end
    
    %% Applying free surface
    txx(2:M+1,2:N+1)=txx(2:M+1,2:N+1)+lambdaplus2mudtbydh*(u(2:M+1,2:N+1)-u(2:M+1,1:N))+lambdadtbydh*(v(2:M+1,2:N+1)-v(1:M,2:N+1));
    tzz(2:M+1,2:N+1)=tzz(2:M+1,2:N+1)+lambdaplus2mudtbydh*(v(2:M+1,2:N+1)-v(1:M,2:N+1))+lambdadtbydh*(u(2:M+1,2:N+1)-u(2:M+1,1:N));
    txz=txz+dtmubydh*(v(1:M+1,2:N+2)-v(1:M+1,1:N+1)+u(2:M+2,1:N+1)-u(1:M+1,1:N+1));

    %% Generate resultant velocity matrix and image plot it
    resvty=u(2:M+1,1:N).^2+v(2:M+1,1:N).^2;
     imagesc(resvty,'parent',uaxis);
    title(strcat('i=',num2str(tind),'maxval=',num2str(max(max(resvty)))));
    axis('equal');
    colormap(gray);
    pause(0.001);
end


