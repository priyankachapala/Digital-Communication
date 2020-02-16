%Hadamard codes
clc;

clear all;
N=32;

H32 = ones(N);
k = N/2;
while k>=1
    
    for i=1:N/k
        for j=1:N/k
            
            if mod(i,2)+mod(j,2)==0
                H32((i-1)*k+1:i*k,(j-1)*k+1:j*k)=~H32((i-1)*k+1:i*k,(j-1)*k+1:j*k);
                
            end
            
        end
        
    end
    k=k/2;
    
end

H32 = 2*H32-1;

figure; subplot(211);

cor([H32(7,:)], [H32(7,:)]); title('Auto-correlation of 32 length Walsh-Hadamard code');

subplot(212);

cor([H32(3,:)],[H32(15,:)]); title('Cross-correlation of 32 length Walsh-Hadamard code');

