function [w, m, ent] = adaptiveEnvelope(signal, fs, varargin)

% [w, m, ent] = adaptiveEnvelope(signal, fs, varargin)
%
% Extracts the RMS EMG envelope using the automatic adaptive procedure
% described in Ranaldi S., De Marchis C., Conforto S. "An automatic,
% adaptive, information-based algorithm for the extraction of the sEMG
% envelope".
%
% [w, m, ent] = adaptiveEnvelope(signal, fs)
%   Extract the envelope from the EMG signal sampled at fs Hz.
%
% [w, m, ent] = adaptiveEnvelope(..., 'mincontrol', val)
%   If val = true sets the minimum of the extracted envelope to 0.
% [w, m, ent] = adaptiveEnvelope(..., 'plotresults', val)
%   If val = true plots the results.
% [w, m, ent] = adaptiveEnvelope(..., 'language', l)
%   If l = 'C' uses the mex-C version.
%   If l = 'MATLAB' uses the MATLAB version.
%
% OUTPUTS:
% 	w: the extracted envelope.
% 	m: the point by point optimal filter lengths.
% 	ent: the estimation entropy at each iteration step.


%% Input control

narginchk(1,7);

minControl=false;
varPlot=false;
language='C';
whitenWindow=max([fs, round(length(signal)/5)]);
if ~isempty(varargin)
    i=1;
    while i<=length(varargin)-1
        ctrl=0;
        switch varargin{i}
            case 'mincontrol'
                if ~islogical(varargin{i+1})
                    error('Invalid argument type');
                end
                minControl=varargin{i+1};
                ctrl=1;
            case 'plotresults'
                if ~islogical(varargin{i+1})
                    error('Invalid argument type');
                end
                varPlot=varargin{i+1};
                ctrl=1;
            case 'language'
                language=varargin{i+1};
                ctrl=1;
            case 'wwin'
                whitenWindow=varargin{i+1};
                ctrl=1;
        end
        if ~ctrl
            error([varargin{i}, ' is not a valid argument']);
        end
        i=i+2;
    end
end


%% Definition of the parameters of the algorithm.

alpha=1;
nu=2;
maxIter=20;

p=(2^(1/2*alpha))*gamma((alpha+1)/(2*alpha))/sqrt(pi); % Normalization Factor
%val=(((sqrt(pi)*gamma(nu+0.5))/(gamma(nu+0.5)))^2 - 1)/((alpha*nu)^2);

%% Initialization.

load('chiTable.mat','chiTable'); % Loads the lookup table for the chi squared distribution.

if size(signal,1)>size(signal,2)
    signal=signal';
end

[signal, signalOld]=conditionEMG(signal,language,whitenWindow);

signal(isnan(signal)) = 0.0001*rand(1); % Provisional

convStep=maxIter*ones(size(signal));
ent=zeros(maxIter,length(signal));
idx=1:length(signal);
idx1=1:length(signal);
idxB=zeros(1,length(signal));
ctrl=0;

if whiteTest(signal)
    warning('The signal is not white. Results can be inaccurate');
end

%% Initialization of the window length.

ll=min([length(signal), 10000]);

% Adaptive initialization.
% [Pxx,F] = periodogram(abs(signal)-nanmean(abs(signal)),[],fs,fs);
% [~, locs]= findpeaks(Pxx);
% if length(locs)>1
%     ll=F(locs(2));
% else
%     ll=F(locs(1));
% end
% ll=1/ll;
% ll=round(ll*1000);

m=ones(size(signal)).*(ll);

%% Test for the whiteness of the signal.

% if ~whiteTest(signal)
%     disp(['Signal is white.']);
% else
%     disp('Signal is not white.');
% end

%% First static estimation.

w=staticEstimationW(signal, m, alpha, p);
[d,d2]=staticEstimationD(signal, m, alpha, p);

m=filterLength(w,d,d2,alpha,nu,idx,m);
mProv(1,:)=m;
[w] = envelopeEstimation(signal,m,alpha,nu,idx,w,p);
wProv(1,:)=w;
dProv(1,:)=d;
d2Prov(1,:)=d2;
count=1;

%% Optimal filter length extraction.

% if language == "C"
%     signal = signal';
% end

switch language
    
    case "C"

        [m,w] = loopFunction(signal,w,d,d2,m,chiTable);
        
    case "MATLAB"
        
        while ctrl==0
            
            
            
            mp=filterLength(w,d,d2,alpha,nu,idx,m);
            [wp] = envelopeEstimation(signal,mp,alpha,nu,idx,w,p);
            [dp,dp2]=derivativesEstimation(signal,mp,alpha,nu,idx,d,d2,p);
            
            
            [ent(count,:)]=estEntropy(signal,m,chiTable);
            
            if count>3
                
                e1=ent(count,:)-ent(count-1,:);
                e2=ent(count-1,:)-ent(count-2,:);
                ii=e2-e1<0; % Convergence test.
                convStep(ii)=min(count-1,convStep(ii));
                idxB=idxB | ii;
                idx=find(~idxB);
                
            end
            
            % Updates
            m(idx)=mp(idx);
            w(idx)=wp(idx);
            d(idx)=dp(idx);
            d2(idx)=dp2(idx);
            
            count=count+1;
            
            % Break conditions.
            if length(idx)<0.05*length(signal)
                ctrl=1;
            end
            
            if count>maxIter
                disp('Maximum number of iterations reached!');
                disp([num2str(100*length(idx)/length(w)),'% of the sample did not converge.']);
                ctrl=1;
            end
            
        end
        
end


%% Envelope extraction.

[w] = envelopeEstimation(signal,m,alpha,nu,idx1,w,p);

if minControl
    if length(find(w<(min(w(20:end-20))+0.05*range(w(20:end-20)))))>0.01*length(w)
        w=w-mean(w(w<(min(w(20:end-20))+0.05*range(w(20:end-20))))); % Offset removal.
    end
end

%convStep(convStep==maxIter)=NaN; % Just for visualization purposes.

%% Plot of the results.

if varPlot
    figure
    a(1)=subplot(2,1,1);
    plot(signalOld,'color',0.5*[1 1 1]); hold on; plot(w,'r','linewidth',3); hold off;
    a(2)=subplot(2,1,2);
    plot(m,'b')
    linkaxes(a,'x');
end