function [der1, der2] = derivativesEstimationMat(signal, m, alpha, nu, idx, der1, der2, p)

% [der1 der2] = derivativesEstimation(signal, m, alpha, nu, idx, der1, der2, w)
%
% Used in the adaptiveEnvelope algorithm

%% Initialization.

% p=(2^(1/2*alpha))*gamma((alpha+1)/(2*alpha))/sqrt(pi); % Normalization factor

est=zeros(size(der1));

%% Derivatives computation.

for i=1:length(idx)
    
    k=idx(i);
    
    semiLen=ceil(m(k)/2);
    
    s=signal(max(1,k-semiLen):min(length(signal),k+semiLen));
    s=abs(s);
    
    a=1:length(s);
    a=a-ceil(length(s)/2);
    
    d=sum(a.^2);
    d2=sum(a.^4);
    
    est(k)=sum((a).*(s.^(1/alpha)));
    est2(k)=sum((a.^2).*(s.^(1/alpha)))./(d2.*p) - (d/(d2*p)).*(sum((1-(a.^2).*(d/d2)).*(s.^(1/alpha))))./(length(s)+(d/d2)*d);
    der1(k)=est(k)./(d*p);
    der2(k)=est2(k)*2;
    
end