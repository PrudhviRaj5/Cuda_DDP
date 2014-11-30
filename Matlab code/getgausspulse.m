function op=getgausspulse(h)
    sigma=(420e-6)/6;
    m=0;
    x=[m-3*sigma:h:m+3*sigma];
    op=exp(-((x-m).^2)/(2*sigma^2));
end