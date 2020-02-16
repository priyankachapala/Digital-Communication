N=10^4;
ip=rand(1,N)>.5;
s=2.*ip-1;
n=1/sqrt(2)*[randn(1,N) + j*randn(1,N)];
eb=[-3:10];
err=[];
for d=1:length(eb)
    y=s+10^(-1*eb(d)/20)*n;
    iphard=real(y)>0;
    err(d)=size(find([ip-iphard]),2);
end
simulatedber=err/N;
theoryber=0.5*erfc(sqrt(10.^(eb/10)));
close all
semilogy(eb,theoryber,'g-');
hold on
semilogy(eb,simulatedber,'b--');
axis([-3 11 10^-5 0.5]);
grid on;
legend('Theoretical','Simulation');
xlabel('Eb/No, dB');
ylabel('bit error rate');
title('Bit Error Prob curve for BPSK modulation');

scatterplot(y);