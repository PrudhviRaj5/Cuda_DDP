clear all;
close all;
clc;
p=2667;
lo=6396;
sh=3103;
mu=sh^2*p;
lambda=lo^2*p-2*mu;
freq=2.25*1e6;
wavelength=lo/freq;
dh=wavelength/30;
dt=dh/(lo*1.5);
length=20;                   %in mm
breadth=20;                  %in mm
xno=floor(length*1e-3/dh);
yno=floor(breadth*1e-3/dh);
xno=xno+16-mod(xno,16);
yno=yno+16-mod(yno,16);
xno=xno-2;
yno=yno-2;
xsource=floor((xno+2)/2);
ysource=floor((yno+2)/2);
time_total=10e-6;
timesteps=ceil(time_total/dt);
 clrmap=load('rainbowcolormap.txt');
X=zeros(xno+2,yno+2);
Y=zeros(xno+2,yno+2);
Z=zeros(xno+1,yno+1);
U=zeros(xno+2,yno+2);
V=zeros(xno+2,yno+2);
x=zeros(xno+2,yno+2);
y=zeros(xno+2,yno+2);
for n=1:timesteps
%      applypml(X,0,20,20,20);
%      applypml(Y,20,20,20,20);
%      applypml(Z,20,20,20,20);
    U(2:xno+1,1:yno+1)=U(2:xno+1,1:yno+1)+(dt/(p*dh))*(X(2:xno+1,2:yno+2)-X(2:xno+1,1:yno+1)+Z(2:xno+1,1:yno+1)-Z(1:xno,1:yno+1));
    V(1:xno+1,2:yno+1)=V(1:xno+1,2:yno+1)+(dt/(p*dh))*(Z(1:xno+1,2:yno+1)-Z(1:xno+1,1:yno)+Y(2:xno+2,2:yno+1)-Y(1:xno+1,2:yno+1));
      if n*dt<=3/freq
            V(ysource,xsource)=(1-cos(2*pi*freq*dt*n/3))*cos(2*pi*freq*n*dt);
      end
   
    X(2:xno+1,2:yno+1)=X(2:xno+1,2:yno+1)+((lambda+2*mu)*(dt/dh)*(U(2:xno+1,2:yno+1)-U(2:xno+1,1:yno)))+((lambda*dt)/dh)*((V(2:xno+1,2:yno+1)-V(1:xno,2:yno+1)));
    Y(2:xno+1,2:yno+1)=Y(2:xno+1,2:yno+1)+((lambda+2*mu)*(dt/dh)*(V(2:xno+1,2:yno+1)-V(1:xno,2:yno+1)))+((lambda*dt)/dh)*((U(2:xno+1,2:yno+1)-U(2:xno+1,1:yno)));
    Z=Z+((mu*dt)/dh)*(V(1:xno+1,2:yno+2)-V(1:xno+1,1:yno+1)+U(2:xno+2,1:yno+1)-U(1:xno+1,1:yno+1));
   x(:,:)=x(:,:)+U(:,:)*dt;
   y(:,:)=y(:,:)+V(:,:)*dt;
   result=U.^2+V.^2;
    Result=sqrt(x.^2+y.^2);   
   
%     subplot(211);
    q = mesh(1:xno+2,1:yno+2,Result);hold on;colormap hsv;%colorbar;
    axis([0 xno+2 0 yno+2 -(1e-7) 1e-7]);
%     subplot(212);scatter(n,Result,5,[1 0 0],'fill');axis([0 timesteps -(1e-7) (1e-7)]);hold on;
    title(strcat('timesteps=',num2str(n),' maxvalue=',num2str(max(max(Result)))));
   % M(n)=getframe;
   if (n~=timesteps)
    pause(0.001);delete(q);
   end
 end
 %   plot(1:n,z);