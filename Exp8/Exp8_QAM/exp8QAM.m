clear all
close all
clc;

snrdB = [0:2:12];
snr = 10.^(snrdB/10);%snr conversion from db to normal scale
%3 types of constellations, rectangle, circular 1 and 2
ber_r=zeros(1,length(snrdB));
ber_c1=zeros(1,length(snrdB));
ber_c2=zeros(1,length(snrdB));

ber_rth=zeros(1,length(snrdB));
ber_c1th=zeros(1,length(snrdB));
ber_c2th=zeros(1,length(snrdB));

N=30000;%3 bits*10000 symbols
iter_max=50;
% constellations
rec=[-3+1i -1+1i 1+1i 3+1i -3-1i -1-1i 1-1i 3-1i];
cir1=[1 -1 1i -1i (1+1i)/sqrt(2) (-1+1i)/sqrt(2) (1-1i)/sqrt(2) (-1-1i)/sqrt(2)];
cir2 = [1+1i 1-1i -1+1i -1-1i 1+sqrt(3) -(1+sqrt(3)) (1+sqrt(3))*1i -(1+sqrt(3))*1i];

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
    symbol=signal(m:m+2);% group 3 bits as a symbol
%Check to what symbol it corresponds and assign to it properly
    if(symbol == [0 0 0])
        const_rec((m-1)/3+1)=-3+1i;%INDEXING PROPERLY, MATLAB from 1
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
%Add Additive white gaussian noise
received_rec=awgn(const_rec, snrdB(i)+4.771, 'measured');
received_cir1=awgn(const_cir1, snrdB(i)+4.771, 'measured');
received_cir2=awgn(const_cir2, snrdB(i)+4.771, 'measured');

for l=1:1:length(received_rec)
    dist_rec = zeros(1, 8);
    dist_cir1 = zeros(1, 8);
    dist_cir2 = zeros(1, 8);
    %calculate distance of every symbol with these 8 fixed symbols 
    for ll=1:1:8
        dist_rec(ll)=abs(received_rec(l)-rec(ll));
        dist_cir1(ll)=abs(received_cir1(l)-cir1(ll));
        dist_cir2(ll)=abs(received_cir2(l)-cir2(ll));
    end
   
    %take min distance out of these values
    [v, minrec] = min(dist_rec);
    [v, mincir1] = min(dist_cir1);
    [v, mincir2] = min(dist_cir2);
    %assign that closest symbol codeword to decoder
    decoded_rec(l, :) = rec_bit(minrec, :);
    decoded_cir1(l, :) = cir1_bit(mincir1, :);
    decoded_cir2(l, :) = cir2_bit(mincir2, :);
end
 
    temp = decoded_rec';
    rec_stream = temp(:);%alligning all the symbols(bits) in a row
    temp = decoded_cir1';   
    cir1_stream = temp(:);
    temp = decoded_cir2';
    cir2_stream = temp(:);
    %Calculate error i.e if decoded different from input
    % and taking average over all the iterations
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
title('Constellation of rectangular QAM');
scatterplot(const_cir1);
xlabel('In-phase');
ylabel('Quadrature');
title('Constellation of Circular QAM - 1');
scatterplot(const_cir2);
xlabel('In-phase');
ylabel('Quadrature');
title('Constellation of Circular QAM - 2');
scatterplot(received_rec);
xlabel('In-phase');
ylabel('Quadrature');
title('Constellation of received rectangular QAM');
scatterplot(received_cir1);
xlabel('In-phase');
ylabel('Quadrature');
title('Constellation of received Circular QAM - 1');
scatterplot(received_cir2);
xlabel('In-phase');
ylabel('Quadrature');
title('Constellation of received Circular QAM - 2');
figure;
semilogy(snrdB, ber_r, 'r', snrdB, ber_c1, 'g', snrdB, ber_c2, 'b');
xlabel('Eb/No in dB');
ylabel('Bit Error Rate');
legend('Rectangular constellation', 'Circular constellation 1', 'Circular constellation 2');
title('Comparison of the three schemes, experimental');
figure;
semilogy(snrdB, ber_rth, 'r', snrdB, ber_c1th, 'g', snrdB, ber_c2th, 'b');
xlabel('Eb/No in dB');
ylabel('Bit Error Rate');
legend('Rectangular constellation', 'Circular constellation 1', 'Circular constellation 2');
title('Comparison of the three schemes, theoritical');
figure;
semilogy(snrdB, ber_r, 'r', snrdB, ber_rth, 'r*');
xlabel('Eb/No in dB');
ylabel('Bit Error Rate');
legend('Simulated', 'Theoretical');
title('Rectangular constellation');
figure;
semilogy(snrdB, ber_r, 'g', snrdB, ber_rth, 'g*');
xlabel('Eb/No in dB');
ylabel('Bit Error Rate');
legend('Simulated', 'Theoretical');
title('Circular constellation 1');
figure;
semilogy(snrdB, ber_r, 'b', snrdB, ber_rth, 'b*');
xlabel('Eb/No in dB');
ylabel('Bit Error Rate');
legend('Simulated', 'Theoretical');
title('Circular constellation 2');
 