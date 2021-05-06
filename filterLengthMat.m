function m = filterLengthMat(w,d,d2,alpha,nu,idx,m)

% m = filterLength(w,d,d2,alpha,nu,idx,m)
%
% Used in the adaptiveEnvelope algorithm

%% Initialization.

aa=zeros(size(w));
bb=zeros(size(w));
num=zeros(size(w));
den=zeros(size(w));

%% Derivatives definition.

C=4*f(alpha,nu);
aa(idx)=-d(idx)/2;
bb(idx)=(1/6) * (d2(idx) + ((alpha*nu - 1)*d(idx).^2)./(4.*w(idx)));

%% Optimal filter length computation.

num(idx)=(C.*(w(idx).^4));
den(idx)=((bb(idx).*w(idx))+(alpha*nu - 1).*(aa(idx).^2))./2;
den(idx)=den(idx).^2;

m(idx)=abs(num(idx)./den(idx));
m(idx)=round(m(idx).^(1/5));

m(m>=10000)=10000;
m(m<1)=1;