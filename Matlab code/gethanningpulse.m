function op=gethanningpulse(freq,dt,tmax)
    for ind=1:tmax
        if ind*dt<=3/freq
            sig(ind)=(1-cos(2*pi*freq*dt*ind/3))*cos(2*pi*freq*ind*dt);
        else
            sig(ind)=123456;
        end
    end
    op=sig;
end