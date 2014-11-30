clc
par=load('simdata');
dely=par.simdata.tdelay;
rcvdelay=max(dely)-dely;
for scanind=0:50
    rcvdata=load(strcat('scan_',num2str(scanind)));
    asc=zeros(1,6000);
    for rcvind=1:16
        delaydasc=rcvdata.rcv(rcvind).u(rcvdelay(rcvind)+1:end);
        asc(1:length(delaydasc))=asc(1:length(delaydasc))+delaydasc;
    end
    bsc(scanind+1,:)=asc;
end
        