function [ePow]=estEntropy(signal,m,chiTable)

% [ePow]=estEntropy(signal,m,chiTable)
%
% Used in the adaptiveEnvelope algorithm

%% Point by point entropy evaluation.

for k=1:length(signal)
    
    semiLen=ceil(m(k)/2);
        
    ns=chiTable(m(k),:);
    
    aa=ns.*log(1./ns);
    aa(ns==0)=0;
    ePow(k)=sum(aa);
        
end



