clear all;
close all;
N = 10^5;
ip = rand(1,N)>0.5;
s=2*ip-1;
f=0:2*pi/100:2*pi-2*pi/100;
c=cos(2*pi*f);
c1 = c/(sqrt(sum(c.^2)));
s1=s.'*c1;
n=(1/sqrt(2))*(randn(N,100)+1i*randn(N,100));
ebn0db=[-3:10];
for ii=1:length(ebn0db)
    y1=s1+10^(-1*ebn0db(ii)/20)*n;
    y=c1*y1';
    iphard = real(y)>0;
    errn(ii)=size(find([ip-iphard]),2);
end
figure
simber=errn/N;
theoryber=0.5*erfc(sqrt(10.^(ebn0db/10)));
close all
semilogy(ebn0db,theoryber,'b');
hold on
semilogy(ebn0db,simber,'r');
axis([-3 11 10^-5 0.5]);
grid on;
legend('theory', 'simulation');
xlabel('Eb/N0,db');
ylabel('bit error rate');
title('BER');
figure
subplot(2,1,1),plot(ip(1:10)),title('input');
subplot(2,1,2),plot(iphard(1:10)),title('demodulated');
figure
scatterplot(iphard);
