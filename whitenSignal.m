function [signalOut,coeff] = whitenSignal(signalIn, predictorOrder, windowLength)

% [signalOut] = whitenSignal(signalIn, predictorOrder, windowLength)
%
% Prewhitens the EMG signal using an autoregressive model.
%
% INPUTS:
% 	signalIn: The EMG signal.
% 	predictorOrder: The order of the AR model.
% 	windowLength: The length of the window to be modelled.
%
% OUTPUTS:
% 	signalOut: The whitened EMG signal.

%% Initialization.

% normFactor=sum(signalIn.^2);
semiLen=ceil(windowLength/2);
    
%% AR filtering.

for k=1:length(signalIn)
    
    s=signalIn(max(1,k-semiLen):min(k+semiLen,length(signalIn)));
    
    coeff=aryule(s,min(predictorOrder,length(s)));
    
    e=filter(-coeff,1,signalIn);
    
    signalOut(k)=e(k);
    
end

% signalOut=(signalOut./sum(signalOut.^2)).*normFactor;