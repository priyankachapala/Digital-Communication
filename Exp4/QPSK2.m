n=50000;
dig_data=randn(1,n)>0.5;
dig_data=2*dig_data-1;
inphase=dig_data(1:2:n-1);
quadrature=dig_data(2:2:n);
data=inphase+j*quadrature;
Eb=1;
SNRdb=0:2:14;
SNR=10.^(SNRdb/10);
N=Eb./SNR;
for count=1:length(SNR)
    N0=Eb/SNR(count);
    awgn=sqrt(N0/2)*randn(1,length(inphase))+j*sqrt(N0/2)*randn(1,length(quadrature));
    y=data+awgn;
    rcvd_inphase=real(y);
    rcvd_quadrature=imag(y);
    err=0;
    err1=0;
    for i=1:length(inphase)
        if((rcvd_inphase(i)>0 && inphase(i)==-1)||(rcvd_inphase(i)<0 && inphase(i)==1))
            err=err+1;
        end 
    end
    
    for i=1:length(quadrature)
        if((rcvd_quadrature(i)>0 && quadrature(i)==-1)||(rcvd_quadrature(i)<0 && quadrature(i)==1))
            err=err+1;
        end 
    end
    BER_sim(count)=err/(length(inphase)+length(quadrature));
    for i=1:length(inphase)
        if((rcvd_inphase(i)>0 && inphase(i)==-1) ||(rcvd_inphase(i)<0 && inphase(i)==1)||(rcvd_quadrature(i)>0 && quadrature(i)==-1)||(rcvd_quadrature(i)<0 && quadrature(i)==1))
            err1=err1 + 1;
        end 
    end
    SER_sim(count)=err1/(length(inphase));
end
BER_ideal=(1/2)*erfc(sqrt(SNR));
SER_ideal=2*BER_ideal;
scatterplot(data);
axis([-2 2 -2 2])
title('Constellation diagram without noise');
scatterplot(y);
axis([-2 2 -2 2])
title(' Constellation disgram with noise');
figure,semilogy(SNRdb, BER_ideal, 'k');
hold on
semilogy(SNRdb, SER_sim,'r*');
legend('Ideal', 'Simulated', 3);
axis([min(SNRdb) max(SNRdb) 10^(-5) 1]);
title('BER');
hold off
 
figure,semilogy(SNRdb, SER_ideal, 'k');
hold on
semilogy(SNRdb, SER_sim,'r*');
legend('Ideal', 'Simulated', 3);
axis([min(SNRdb) max(SNRdb) 10^(-5) 1]);
title('SER');
hold off
 
h=rcosdesign(0.6,10,4);
fvtool(h,'analysis','impulse');
x1=upsample(data,4);
re_x=conv(real(x1),h);
im_x=conv(imag(x1),h);
x=re_x + j*im_x;
eyediagram(x,100);
rcvd_noise=rcvd_inphase+ j*rcvd_quadrature;
z1=upsample(rcvd_noise, 4);
z=conv(z1,h);
eyediagram(z,100);
 %hmf=phased.MatchedFilter('Coefficients',getMatchedFilter(h));
 %rcvd_m=stem(hmf,2);
 for count=1:length(SNR)
     N0=Eb/SNR(count);
     awgn=sqrt(N0/2)*randn(1,length(x))+j*sqrt(N0/2)*randn(1,length(x));
     y2=x + awgn;
     re_mfo=conv(real(y2),h);
     im_mfo=conv(imag(y2),h);
     mfo=re_mfo + j*im_mfo;
     mfo_down=downsample(mfo(41:length(mfo)),4)
     rcvd_inphase=real(mfo_down);
     rcvd_quadrature=imag(mfo_down);
     err=0;
     err1=0;
     for i=1:length(data)
         if(((rcvd_inphase(i)>0) && (real(data(i))==-1)) ||((rcvd_inphase(i)<0) &&(real(data(i))==1)))
             err=err+1;
         end
     end
     for i=1:length(data)
         if((rcvd_quadrature(i)>0 && imag(data(i))== -1) || (rcvd_quadrature(i) <0 && imag(data(i))==1))
             err= err+1;
         end
     end
     BER_sim(count)= err/(length(data) + length(data));
     for i=1:length(data)
         if ((rcvd_quadrature(i)>0 && imag(data(i))== -1) || (rcvd_quadrature(i)<0 && imag(data(i))== 1))
             err=err+1;
         end
     end
     BER_sim(count)= err/(length(data) + length(data));
     for i=1:length(data)
         if((rcvd_inphase(i)>0 && real(data(i))== -1) || (rcvd_inphase(i)<0 && real(data(i))== 1) || (rcvd_quadrature(i)>0 && imag(data(i))== -1) || (rcvd_quadrature(i)<0 && imag(data(i))== 1))
             err1= err1+1;
         end
     end
     SER_sim(count)= err1/(length(data));
 end
 BER_ideal= (1/2)*erfc(sqrt(SNR));
 SER_ideal= 2*BER_ideal;
 
 scatterplot(data);
 axis([-2 2 -2 2])
 title(' Constellation diagram without noise - Pulse shaping ');
 scatterplot(y);
 axis([-2 2 -2 2])
 title(' Constellation diagram with noise - Pulse shaping');
 figure, semilogy(SNRdb, BER_ideal, 'k');
 hold on
 semilogy(SNRdb, BER_sim, 'b*');
 legend(' Ideal ', 'Simulated', 3);
 axis([min(SNRdb) max(SNRdb) 10^(-5) 1]);
 title(' BER - Pulse shaping');
 hold off
 
 figure, semilogy(SNRdb, SER_ideal, 'k');
 hold on
 semilogy(SNRdb, SER_sim, 'r*');
 legend(' Ideal ', 'Simulated', 3);
 axis([min(SNRdb) max(SNRdb) 10^(-5) 1]);
 title(' SER - Pulse shaping');
 hold off
 SNRdb= 0:2:14;
 SNR= 10.^(SNRdb/10);
 
 
 
    
  

