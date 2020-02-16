clear all; 
ctr=0; 
n=15; 
k=7; 
t=2; 
N=k*10000; 
ebyn0db = -10 : 2 : 10; 
ebyn0 = 10.^(ebyn0db/10); 
ip=rand(1,N)>0.5; 
si=2*ip-1; 
g=[1 0 0 0 1 0 1 1 1 0 0 0 0 0 0]; 
gen=[]; 
for i=0:k-1 
gg= circshift(g',i); 
gen=[gen gg]; 
end 
gen=gen'; 
a=0:15; 
b=gf(a,4); 
alpha=b(3); 
A=[]; 
for i= 0: 15 
aa= alpha^i; 
A=[A aa]; 
end 
s=reshape(ip,N/k,k); 
coded = mod(s*gen,2); 
code = reshape(coded',1,(N/k)*n); 
BER = zeros(1,length(ebyn0db)); 
nn = randn(1,N) + 1i*randn(1,N); 
nn1 = randn(N/k,n) + 1i*randn(N/k,n); 
for ii=1: length(ebyn0db) 
rxd1 = si + 10^(-ebyn0db(ii)/20)*nn; 
rxd = real(rxd1)>0; 
nerr(ii)=size(find(ip-rxd),2); 
coded_n = (2*coded-1) + 10^(-ebyn0db(ii)/20)*nn1; 
coded_n = real(coded_n) >0; 
rcv=[]; 
for i=1:N/k 
pos=find(coded_n(i,:)); 
for m=1:length(pos); 
if m==1 
s1= A(pos(m)); 
s3= A(pos(m))^3; 
else
s1=s1+A(pos(m)); 
s3=s3+(A(pos(m)))^3; 
end 
end 
s1d=double(s1.x);
s11=de2bi(s1d,2*t); 
s3d=double(s3.x); 
if (s1==b(1)) && (s3==b(1)) 
rcv=[rcv coded_n(i,:)]; 
elseif (s1 ~= b(1)) && (s3 == s1^3) 
e=s11; 
im = find(A(1:15) == s1); 
rr=xor(coded_n(i,:),[zeros(1,im-1) 1 zeros(1,n-im)]); 
rcv=[rcv rr]; 
elseif (s1 ~= b(1)) && (s3 ~= s1^3) 
pr=(s3+s1^3)/s1; er=[]; 
for f=1:length(b) 
if (b(f))^2 + s1*b(f) + pr == b(1) 
er = [er b(f)]; 
end 
end 
if isempty(er) 
rcv=[rcv coded_n(i,:)]; 
ctr=ctr+1; 
else e1 = find(A(1:15) == er(1)); 
e2 = find(A(1:15) == er(2)); 
if(e1>e2) 
xx= [zeros(1,e2-1) 1 zeros(1,e1-1-e2) 1 zeros(1,n-e1)]; 
elseif(e1<e2) 
xx= [zeros(1,e1-1) 1 zeros(1,e2-1-e1) 1 zeros(1,n-e2)]; 
else xx= [zeros(1,e1-1) 1 zeros(1,n-e1-1)]; 
end 
rr= xor(coded_n(i,:),xx); 
rcv =[rcv rr]; 
end 
else rcv=[rcv coded_n(i,:)]; 
end 
end 
rcv=reshape(rcv, n, N/k); 
y= length(find(xor(rcv,coded'))); 
BER(ii)= y/(n*(N/k)); 
y1= length(find(xor(rxd,ip))); 
thber_uncoded(ii) = y1/N; 
end 
uncodber=nerr/N; 
thber = 0.5*erfc(sqrt(10.^(ebyn0db/10))); 
figure; 
semilogy(ebyn0db,BER,'md-','Linewidth',3); 
hold on; 
semilogy(ebyn0db,uncodber,'b.-','Linewidth',3); 
hold on; 
semilogy(ebyn0db,thber,'r.-','Linewidth',3); 
hold on; 
legend('Coded AWGN','Uncoded AWGN','THEORETICAL BER'); 
title('BER PLOTS OF CODED AND UNCODED BPSK');
