function Sobs = sm2sobs(Sm, Jv, k)
% SM2SOBS Evaluates the observed sieving coefficients
% from the values of intrinsic sieving coefficients.
%   SOBS = SM2SOBS(SM, JV, K) evaluates the observed 
%   sieving coefficient (SOBS) from the values of the
%   the intrinsic coefficient (SM), filtrate flux (JV) 
%   and mass transfer coefficient (K). The film model
%   is used in the calculations. JV and K must be in 
%   the same unit System.

Sobs = Sm ./ (Sm + (1 - Sm) .* exp(-Jv ./ k));