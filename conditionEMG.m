function [out,outOld]=conditionEMG(signal, language, whitenWindow)

% [out]=conditionEMG(signal)
%
% Conditioning block for an EMG signal.
% Consists of a [20, 450] Hz band pass filter, a 50 Hz notch filter
% and a prewhitening AR filter.
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
        
        out = [];
        k=0;
        
        while k+whitenWindow<=length(signal)
            
            tmp_out = whiteningSignal(signal(k+1:k+whitenWindow), 13, 150);
            out = [out; tmp_out];
            
            k = k + whitenWindow;
            
        end
        
        if k < length(signal)
            out = [out; whiteningSignal(signal(k+1:end), 13, 150)];
        end
        
        out=((out./max(out)).*normFactor)';
        
    case 'MATLAB'
        
        out = [];
        k=0;
        
        while k+whitenWindow<=length(signal)
            
            out = [out, whitenSignal(signal(k+1:k+whitenWindow), 13, 150)];
            
            k = k + whitenWindow;
            
        end
        
        if k < length(signal)
            out = [out, whitenSignal(signal(k+1:end), 13, 150)];
        end
        
        out=((out./max(out)).*normFactor);
        
end