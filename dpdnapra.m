function D = dpdnapra(T, eta, nbp)
% DPDNAPRA evaluates the diffusion coefficient of pDNA by the
% proposed correlation of Prazeres [1].
%   D = DPDNAPRA(T, ETA, NBP) evaluates the diffusion coefficient
%   of pDNA molecules (D) with NBP base pairs, at a temperature
%   of T degrees Kelvin and in a solution with viscosity
%   ETA ([Pa.s]).
%   SI units must be used.
%     
%   References
% 
%   [1] D. M. F. Prazeres, Prediction of diffusion coefficients
%       of plasmids, Biotechnology and Bioengineering 99 (2008)
%       1040-1044.
 
D = 3.31e-15 * T .* (nbp) .^ (-2 / 3) ./ eta;