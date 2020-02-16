clear all;
close all;
N=10^6;
ip = rand(1,N)>0.5;
s=2*ip-1;
n=1/sqrt(2)*[randn(1,N)+j*randn(1,N)];
ebnodb=[-3:10];
errn=[];
for ii=1:length(ebnodb)
    y=s+10^(-1*ebnodb(ii)/20)*n;
    iphard=real(y)>0;
    errn(ii)=size(find([ip-iphard]),2);
end
simber=errn/N;
theoryber=0.5*erfc(sqrt(10.^(ebnodb/10)));
close all
semilogy(ebnodb,theoryber,'b');
hold on
semilogy(ebnodb,simber,'r');
axis([-3 11 10^-5 0.5]);
grid on;
legend('theory', 'simulation');
xlabel('Eb/N0,db');
ylabel('bit error rate');
title('BER');
figure

scatterplot(s);
