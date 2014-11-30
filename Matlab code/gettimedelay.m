function op=gettimedelay(ip)
%Ref:"Beam focusing behavior of linear phased arrays" L. Azar, Y. Shi,
%S.-C. Wooh*
    n=1:ip.numofelements;
    N_bar=(ip.numofelements-1)/2;
    d=ip.pitch;
    F=ip.focaldepth;
    theta_s=ip.angle*pi/180;    %Angle with normal. Towards left is positive
    c=ip.velocity;

 
    A=(1+(N_bar*d/F)^2 +(2*N_bar*d*sin(theta_s)/F))^(1/2);
    B=(1+((n-N_bar-1)*d/F).^2 - 2.*(n-N_bar-1).*d*sin(theta_s)/F).^(1/2);
    op=(F/c)*(A-B);
          
          
%         time=round(((F/c)*(A-B1))*1.6/dt); % to take care of the rounding
%         operation

%     numofelements=ip.numofelements;
%     nbar=(numofelements-1)/2;
%     d=ip.pitch;
%     F=ip.focaldepth;
%     theta=deg2rad(ip.angle);    %Angle with normal. Towards left is positive
%     c=ip.velocity;
%     for n=0:numofelements-1
%         tt1=
%         term1=1+(nbar*d/F)^2 + (2*nbar*d*sin(theta))/F;
%         term2=1+((n-nbar)*d/F)^2 - (2*(n-nbar)*d*sin(theta))/F;
%         res(n+1)=(F/c)*((term1^0.5)-(term2^0.5));
%     end
%     op=res;
end
        