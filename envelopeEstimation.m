function w=envelopeEstimation(signal, m, alpha, nu, idx, w, p)

% w=envelopeEstimation(signal, m, alpha, nu, idx, w)
%
% Used in the adaptiveEnvelope algorithm

%% Initialization.

% p=(2^(1/2*alpha))*gamma((alpha+1)/(2*alpha))/sqrt(pi); % Normalization factor

est=zeros(size(w));

%% Envelope computation.

for i=1:length(idx)
    
    k=idx(i);
    
    semiLen=ceil(m(k)/2);
    
    s=signal(max(1,k-semiLen):min(length(signal),k+semiLen));
    
    est(k)=sum((abs(s)).^nu);
    est(k)=est(k)/length(s);
    w(k)=(est(k)./p).^(1/(alpha*nu));
    
end

