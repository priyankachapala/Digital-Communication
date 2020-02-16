 clear all;
close all;
clc;
 
snrdB = [0:2:18];
snr = 10.^(snrdB/10);
ber_r=zeros(1, length(snrdB));
ber_c1=zeros(1, length(snrdB));
ber_c2=zeros(1, length(snrdB));
 
ber_rth = zeros(1, length(snrdB));
ber_c1th=zeros(1, length(snrdB));
ber_c2th=zeros(1, length(snrdB));
 
N=30000; 
iter_max = 50;
 
rec = [-3+j -1+j 1+j 3+j -3-j -1-j 1-j 3-j];
cir1=[1 -1 j -j (1+j)/sqrt(2) (-1+j)/sqrt(2) (1-j)/sqrt(2) (-1-j)/sqrt(2)];
cir2 = [1+j 1-j -1+j -1-j 1+sqrt(3) -(1+sqrt(3)) (1+sqrt(3))*j -(1+sqrt(3))*j];
 
rec_bit = [0 0 0; 0 1 0; 1 0 0; 1 1 0; 0 0 1; 0 1 1; 1 0 1; 1 1 1];
cir1_bit=[1 1 1; 0 0 1; 0 1 0; 1 0 0; 1 1 0; 0 1 1; 1 0 1; 0 0 0];
cir2_bit=[0 0 0; 0 1 1; 0 0 1; 0 1 0; 1 0 0; 1 1 0; 1 0 1; 1 1 1];
 
for i=1:1:length(snrdB)
    
error = zeros(1, iter_max);
for k=1:1:iter_max
    
    signal = round(rand(1,N));
    const_rec = zeros(1, N/3);
    const_cir1 = zeros(1, N/3);
    const_cir2 = zeros(1, N/3);
    
for m=1:3:length(signal)
    symbol=signal(m:m+2);
    if(symbol == [0 0 0])
        const_rec((m-1)/3+1)=-3+1i;
        const_cir1((m-1)/3+1)=(-1-1i)/sqrt(2);
        const_cir2((m-1)/3+1)=1+1i;
    elseif(symbol==[0 0 1])
        const_rec((m-1)/3+1)=-3-1i;
        const_cir1((m-1)/3+1)=-1;
        const_cir2((m-1)/3+1)=-1+1i;
    elseif(symbol==[0 1 0])
        const_rec((m-1)/3+1)=-1+1i;
        const_cir1((m-1)/3+1)=1i;
        const_cir2((m-1)/3+1)=-1-1i;
    elseif(symbol==[0 1 1])
        const_rec((m-1)/3+1)=-1-1i;
        const_cir1((m-1)/3+1)=(-1+1i)/sqrt(2);
        const_cir2((m-1)/3+1)=1-1i;
    elseif(symbol==[1 0 0])
        const_rec((m-1)/3+1)=1+1i;
        const_cir1((m-1)/3+1)=-1i;
        const_cir2((m-1)/3+1)=1+sqrt(3);
    elseif(symbol==[1 0 1])
        const_rec((m-1)/3+1)=1-1i;
        const_cir1((m-1)/3+1)=(1-1i)/(sqrt(2));
        const_cir2((m-1)/3+1)=(1+sqrt(3))*1i;
    elseif(symbol==[1 1 0])
        const_rec((m-1)/3+1)=3+1i;
        const_cir1((m-1)/3+1)=(1+1i)/sqrt(2);
        const_cir2((m-1)/3+1)=-(1+sqrt(3));
    elseif(symbol==[1 1 1])
        const_rec((m-1)/3+1)=3-1i;
        const_cir1((m-1)/3+1)=1;
        const_cir2((m-1)/3+1)=-(1+sqrt(3))*1i;
    end
end
 
received_rec=awgn(const_rec, snrdB(i)+4.771, 'measured');
received_cir1=awgn(const_cir1, snrdB(i)+4.771, 'measured');
received_cir2=awgn(const_cir2, snrdB(i)+4.771, 'measured');
 
for l=1:1:length(received_rec)
    dist_rec = zeros(1, 8);
    dist_cir1 = zeros(1, 8);
    dist_cir2 = zeros(1, 8);
    
    for ll=1:1:8
        dist_rec(ll)=abs(received_rec(l)-rec(ll));
        dist_cir1(ll)=abs(received_cir1(l)-cir1(ll));
        dist_cir2(ll)=abs(received_cir2(l)-cir2(ll));
    end
    
    [v, minrec] = min(dist_rec);
    [v, mincir1] = min(dist_cir1);
    [v, mincir2] = min(dist_cir2);
    
    decoded_rec(l, :) = rec_bit(minrec, :);
    decoded_cir1(l, :) = cir1_bit(mincir1, :);
    decoded_cir2(l, :) = cir2_bit(mincir2, :);
end
 
    temp = decoded_rec';
    rec_stream = temp(:);
    temp = decoded_cir1';   
    cir1_stream = temp(:);
    temp = decoded_cir2';
    cir2_stream = temp(:);
    
    error_rec = abs(rec_stream'-signal);
    errorno_rec(k) = sum(error_rec);
    ber_r(i) = ber_r(i) + (errorno_rec(k)/iter_max);
        
    error_cir1 = abs(cir1_stream'-signal);
    errorno_cir1(k) = sum(error_cir1);
    ber_c1(i) = ber_c1(i) + (errorno_cir1(k)/iter_max);
        
    error_cir2 = abs(cir2_stream'-signal);
    errorno_cir2(k) = sum(error_cir2);
    ber_c2(i) = ber_c2(i) + (errorno_cir2(k)/iter_max);
end
 
end
 
ber_r = ber_r/N;
ber_c1 = ber_c1/N;
ber_c2 = ber_c2/N;
 
ber_rth = 2*erfc(sqrt((9/14)*snr));
ber_c1th = erfc(sqrt(3*snr/14))*sin(pi/8);
ber_c2th= 3.5*erfc(sqrt(snr));
scatterplot(const_rec);
xlabel('In-phase');
ylabel('Quadrature');
title('Constellation diagram of rectangular 8-QAM(Theoretical)');
scatterplot(const_cir1);
xlabel('In-phase');
ylabel('Quadrature');
title('Constellation diagram of Circular-1  8-QAM (Theoretical)');
scatterplot(const_cir2);
xlabel('In-phase');
ylabel('Quadrature');
title('Constellation diagram of Circular-2 8-QAM (Theoretical)');
scatterplot(received_rec);
xlabel('In-phase');
ylabel('Quadrature');
title('Constellation diagram of rectangular 8-QAM at decoder');
scatterplot(received_cir1);
xlabel('In-phase');
ylabel('Quadrature');
title('Constellation diagram of received Circular-1 8-QAM at decoder');
scatterplot(received_cir2);
xlabel('In-phase');
ylabel('Quadrature');
title('Constellation diagram of received Circular-2 8-QAM at decoder');
 
figure;
semilogy(snrdB, ber_r, 'k*-', snrdB, ber_c1, 'b-', snrdB, ber_c2, 'c+-');
xlabel('---->  Eb/No(dB)');
ylabel('---->  Bit Error Rate');
legend('Rectangular constellation', 'Circular constellation 1', 'Circular constellation 2');
title('BER plots - Practically observed');
 
figure;
semilogy(snrdB, ber_rth, 'k*-', snrdB, ber_c1th, 'b-', snrdB, ber_c2th, 'c+-');
xlabel('---->  Eb/No(dB)');
ylabel('---->  Bit Error Rate');
legend('Rectangular constellation', 'Circular constellation 1', 'Circular constellation 2');
title('BER plots - Theoretically calculated');


