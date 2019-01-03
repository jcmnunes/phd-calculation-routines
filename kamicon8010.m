function k = kamicon8010(D, w, eta)
% KAMICON8010  evaluates the mass transfer coefficient
% for the Amicon 8010 stirred cell.
%   K = KAMICON8010(D, W, ETA) returns the value of the 
%   mass transfer coefficient (K, [m/s]) from the 
%   input values of the diffusion coefficient (D, [m2/s]),
%   the stirring speed (W [rad/s]) and the solution 
%   viscosity ([ETA, [Pa.s]).
%   The Opong and Zydney's equation is used [1].
%   SI units must be used.
%   The function is vectorized. If D, W, and ETA are
%   vectors, they must be the same size.
%
%   References
%
%   [1] W. S. Opong, A. L. Zydney, Diffusive and convective 
%       protein transport through asymmetric membranes, AIChE
%       Journal 37 (1991) 1497-1510.

a     = 0.23;
b     = 0.567;
c     = 0.33;
ro    = 1000;
rcell = 12.5e-3;

k = D ./ rcell * a .* (w * rcell ^ 2 * ro ./ eta) .^ b .* ...
    (eta ./ (ro * D)) .^ c;