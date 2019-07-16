function [h, p_value, R] = whiteTest(x, maxlag, alpha)

% Tests against the null hyptothesis of whiteness 
% see http://dsp.stackexchange.com/questions/7678/determining-the-whiteness-of-noise for details
%
% demo: 
% % white:
% >> [h,p,R]=white_test(((filter([1], [1], rand(1e3,1)))))
% >> h = 0, p =1, R = 455
% % non white:
% >> [h,p,R]=white_test(((filter([1 .3], [.4 0.3], rand(1e3,1)))))
% >> h = 1, p = 0, R = 2e3

N = length(x);
x = x-mean(x);

if ~exist('m','var')
    maxlag = N-1;
end
if ~exist('alpha','var')
    alpha = 0.05;
end

[r, lag] = xcorr(x, maxlag, 'biased'); 
R = N/r(lag==0)^2*sum(r(lag > 0).^2);
p_value = 1-chi2cdf(R, maxlag);
T = chi2inv(1-alpha, maxlag);
h =  R > T; 