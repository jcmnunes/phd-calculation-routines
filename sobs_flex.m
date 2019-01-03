function [Sobs] = sobs_flex (data)
% SOBS_HS  Calculations of the observed sieving
% coefficient for solutes modeled as FJC,
% in Amicon 8010 filtration cell.
%   [SOBS] = SOBS_HS (DATA) returns the observed 
%   sieving coefficient (SOBS), calculated by the
%   hindered transport model applied to flexible
%   macrosolutes [1]. The model applies to both
%   neutral solutes and membranes.
%   DATA must be a matrix with the following entries
%   (without column headers):
%
%      Jv   w    rp   rs   T   
%      -----------------------
%      num  num  num  num  num  
%      num  num  num  num  num
%      ...
%
%      Jv - filtration flux [m/s] 
%      w  - stirring speed [rad/s]
%      rp - membrane pore radius [m]
%      rs - solute stokes radius [m]
%      T  - temperature [K]
%
%   SI units must be used
%
%   For more information on the hindered transport
%   model, applied to flexible macrosolutes, see [1].
%
%   References
%   
%   [1] Morao et al, Journal of Membrane Science, 336
%       (2009) 61-70
%   
%   see also fiflex, kamicon8010, sm2sobs


Jv = data(:, 1);
w  = data(:, 2);
rp = data(:, 3);
rs = data(:, 4);
T  = data(:, 5);

eta     = visc(T);
rg      = 1.505 * rs;
lambdag = rg ./ rp;

D  = stokes_einstein(rs, eta, T);
Sm = fiflex(lambdag);
k  = kamicon8010(D, w, eta);

Sobs = sm2sobs(Sm, Jv, k);