function [mant, expnt] = mantexpnt(arg)
% MANTEXPNT  returns the mantissa and exponent of a 
% real base 10 argument.
%   [MANT, EXPNT] = MANTEXPNT(ARG) returns the mantissa
%   (MANT) and exponent (EXPNT) of a real base 10 input
%   argument ARG. 
%   ARG must be scalar.

% Author: Vioarr
% Data: 2007

sgn   = sign(arg);
expnt = fix(log10(abs(arg)));
mant  = sgn * 10 ^ (log10(abs(arg)) - expnt);
if abs(mant) < 1
    mant  = mant * 10;
    expnt = expnt - 1;
end