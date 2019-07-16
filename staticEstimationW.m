function [est] = staticEstimationW(signal,m,alpha,p)

% [est] = staticEstimationW(signal,m,alpha,nu)
%
% Used in the adaptiveEnvelope algorithm

%% Estimation.

for k=1:length(signal)
    
    semiLen=ceil(m(k)/2); 
    
    s=signal(max(1,k-semiLen):min(length(signal),k+semiLen));
    
    est(k)=sum(abs(s).^(1/alpha))/length(s);
    
end

%% Normalization

est=est./p;