function [out,outOld]=conditionEMG(signal,language)

% [out]=conditionEMG(signal)
%
% Conditioning block for an EMG signal.
%
% INPUTS:
% 	signal: The EMG signal to be conditioned.
% 	language: Which version to be used.
%         'C': mex-C version.
%         'MATLAB': MATLAB version.
% OUTPUTS:
% 	out: The conditioned EMG signal.

outOld=signal;

normFactor=max(signal);

%% Prewhitening.

switch language
    
    case 'C'
        
        signal = signal';
        
        [out] = whiteningSignal(signal, 13, 150);
        
        out=((out./max(out)).*normFactor)';
        
    case 'MATLAB'
        
        [out] = whitenSignal(signal, 13, 150);
        
        out=((out./max(out)).*normFactor);
        
end