%main code 
clc;
close all;
m =  5 ;
G = [ 1 0 0 1 0 1 0 1 1 1 1 0 1 0];
    
    C1 = pnsq(G(1:7),m);
    
    C2 = pnsq(G(8:14),m); figure;
    
    subplot(211)
    
    y1=cor(C1,C1);
    stem(y1);
    title(['Auto-correlation of ',num2str(2^m-1),'PN -sequence']); subplot(212)
    
    y2=cor(C1,C2);
    stem(y2);
    title(['Cross-correlation of ',num2str(2^m-1),'PN -sequence']);
    
