function argout = stokes_einstein(argin, eta, T)
% STOKES_EINSTEIN  implementation of the
% Sotkes-Einstein equation.
%   ARGOUT = STOKES_EINSTEIN(ARGIN, ETA, T)
%   If ARGIN is the Stokes radius, ARGOUT is the
%   diffusion coefficient. If ARGIN is the
%   diffusion coefficient, ARGOUT is the Stokes
%   radius. ETA is the viscosity of the solution
%   [Pa.s], T is the absolute temperature [K].
%
%   SI units.

% Boltzmann's constant
kB = 1.3806488e-23;
argout = kB * T ./ (6 * pi * eta .*argin);