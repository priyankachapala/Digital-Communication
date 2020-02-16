%Correlation

function [corr] = cor(C1,C2) 
C1 = [C1 C1];

corr = [];

N = length(C2); 
for i=1:N

corr = [corr sum( C1(i:N+i-1).*C2 )/N];

end
stem(corr);
end

