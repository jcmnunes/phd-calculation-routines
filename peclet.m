function Pe = peclet(rp, Kd, Kc, Jv, D, Lp)
% PECLET  Compute Pe number 
%   PE = PECLET(RP, KD, KC, JV, D, LP) evalutes
%   the Peclet number (PE). Input arguments are:
%       * RP - Membrane pore radius [m]
%       * KD - Hindrance diffusion coefficient
%       * KC - Hindrance convection coefficient
%       * JV - Filtration flux [m/s]
%       * D  - Diffusion coefficient
%       * LP - Hydraulic permeability [m]
%   Use consistent units (Peclet number must be
%   adimensional).

Pe = Kc .* rp.^2 .* Jv ./ (8 .* Kd .* D .* Lp);