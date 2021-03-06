fpass=1500;
fstop=1000;
fs=8000;
d1=0.9;
d2=50;
w=(fpass+fstop)*2*pi/(2*fs);
tb=(fpass-fstop)*2*pi/fs;
N=(8*pi)/tb;
h=zeros(1,N);
for n=1:N
    h(n)=-sin(w*(n-ceil(N/2)))/(pi*(n-ceil(N/2)));
end
h(ceil(N/2))=1-(w/pi);
figure
stem(h)
xlabel('n');
ylabel('h(n)');
title('HPF impulse response without windowing');
fvtool(h);
hold on;
title('HPF without windowing');
for x=1:N
    w(x)=0.54-0.46*(1-cos((2*pi*x/(N-1))));
end
figure
stem(w);
title('Hamming window');
xlabel('n');
ylabel('w(n)');
hw=(h.*w);
figure
stem(hw)
xlabel('n');
ylabel('h(n)');
title('Impulse response after windowing');
fvtool(hw);
xlabel('n');
ylabel('h(n)');
title('HPF after windowing');