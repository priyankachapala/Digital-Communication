f=10;
fs=100;
n=1001;
ts=1/fs;
PhasePoints=2*pi*(1:n-1)*ts;
carrier=sin(f*(PhasePoints));
msg=[1 0 1 1 0 1 0 0 1 0];
for j=1:10
    if (msg(j)==0)
        msgp(j)=-1;
    elseif(msg(j)==1)
        msgp(j)=1;
    end
end
for i=1:1000
    extmsgp(i)=msgp(ceil(i/100));
     extmsg(i)=msg(ceil(i/100));
end
for i=1:1000
    bpsk(i)=carrier(i)*extmsgp(i);
end
for m=1:10
    sum=0;
    for n=1:100
        sum=sum+bpsk(100*(m-1)+n);
        demod(m)=sum;
    end
end
figure
subplot(2,1,1),plot(bpsk),title('input');
subplot(2,1,2),plot(demod),title('demodulated');

    
figure
subplot(3,1,1),plot((1:1000),extmsg),title('message'),xlabel('time(sec)'),ylabel('msg');
subplot(3,1,2),plot((1:1000),carrier),title('carrier'),
xlabel('time(sec)'),ylabel('carrier');
subplot(3,1,3),plot((1:1000),bpsk),title('bandpass bpsk'),
xlabel('time(sec)'),ylabel('msg(n)');
