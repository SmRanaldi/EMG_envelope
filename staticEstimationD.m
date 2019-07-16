function [der1 der2] = staticEstimationD(signal, m, alpha, p)

% [der1 der2] = staticEstimationD(signal, m, alpha, nu)
%
% Used in the adaptiveEnvelope algorithm

% p=(2^(1/2*alpha))*gamma((alpha+1)/(2*alpha))/sqrt(pi); % Normalization factor

%% Estimation.

for k=1:length(signal)
    
    semiLen=ceil(m(k)/2);
    
    s=signal(max(1,k-semiLen):min(length(signal),k+semiLen));
    s=abs(s);
    
    a=1:length(s);
    a=a-ceil(length(s)/2);
    
    d(k)=sum(a.^2);
    d2(k)=sum(a.^4);
    
    est(k)=sum((a).*(s.^(1/alpha)));
    est2(k)=sum((a.^2).*(s.^(1/alpha)));
    
end

%% Normalization

der1=est./(d*p);
der2=est2./(d2*p);