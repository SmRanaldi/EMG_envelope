function [val] = f(alpha,nu)

% [val] = f(alpha,nu)
%
% Used in the adaptiveEnvelope function.

val=(((sqrt(pi)*gamma(nu+0.5))/(gamma(nu+0.5)))^2 - 1)/((alpha*nu)^2);