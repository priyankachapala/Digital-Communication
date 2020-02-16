%Gold codes 
clc;

clear ;
m = 5;

G = [1 0 0 1 0 1 0 1 1 1 1 0 1 0];

    
    C11 = (pnsq(G(1:7),m)==1);
    
    C12 = [];                       
    k = 3;
    
    V11 = [C11 C11 C11 C11 C11];
    for j=1:(2^m-1)
        
        C12 = [C12 V11(k)]; 
        k = k+5;
        
    end
    
    X1 = 2*xor(C11,circshift(C12',1)')-1;
    
    X2 = 2*xor(C11,circshift(C12',2)')-1; figure;
    
    
    subplot(211);
    
    cor(X1,X1);
    
    title(['Auto-correlation of ',num2str(2^m-1),' length goldcode']); subplot(212);
    
    cor(X1,X2);
    
    title(['Cross-correlation of ',num2str(2^m-1),' length goldcode']);
