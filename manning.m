function kuhn = manning(I, z)
% MANNING  evaluates the Kuhn distance of DNA molecules
% as a function of the ionic strength of the solution.
%   KUHN = MANNING(I,Z) evaluates the Kuhn distance 
%   (KUHN, [m]) of pDNA molecules as a function of 
%   the ionic strength (I, [mol/m3]) and valence of salt
%   cation (Z). The equation of Manning is used [1].
%   SI units must be used.
%   The function is vectorized. If I, and Z are
%   vectors, they must be the same size.
%
%   References
%
%   [1] G. S. Manning, Biophysical Journal, 91 (2006) 3607-3616.

anull = 7.5e-9;
Rdna  = 1e-9;
lB    = 0.71e-9;
chi   = 4.2;
b     = 0.17e-9;
NA    = 6.02214129e23;
k     = ((8 * pi) * NA * lB * I) .^ 0.5;

kuhn = 2 * (pi * anull / 2) ^ (2 / 3) * Rdna ^ (4 / 3) ./ ...
       (z .^ 2 * lB) .* ((2 * z * chi - 1) .* k * b .* ...
       exp(-k * b) ./ (1 - exp(-k * b)) - 1 - ...
       log(1 - exp(-k * b)));