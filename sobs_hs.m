function [Sobs] = sobs_hs (data)
% SOBS_HS  Calculations of the observed sieving
% coefficient for solutes modeled as hard-spheres,
% in Amicon 8010 filtration cell.
%   [SOBS] = SOBS_HS (DATA) returns the observed 
%   sieving coefficient (SOBS), calculated by the
%   hindered transport model [1], for neutral 
%   membranes and neutral solutes (hard-spheres).
%   DATA must be a matrix with the following entries
%   (without column headers):
%
%      Jv   w    rp   rs   T    Lp
%      -----------------------------
%      num  num  num  num  num  num  
%      num  num  num  num  num  num
%      ...
%
%      Jv - filtration flux [m/s] 
%      w  - stirring speed [rad/s]
%      rp - membrane pore radius [m]
%      rs - solute stokes radius [m]
%      T  - temperature [K]
%      Lp - hydraulic permeability [m]
%
%   SI units must be used
%
%   For more information on the hindered transport
%   model see [1,2].
%
%   References
%
%   [1] W. M. Deen, AIChE Journal, 33 (1987) 1409-1425
%   [2] P. Dechadilok, W. M. Deen, Industrial & Engineering
%       Chemistry Research, 45 (2006) 6953-6959
%   
%   see also kc_dech, kd_dech, sm2sobs, sm_hindered


Jv = data(:, 1);
w  = data(:, 2);
rp = data(:, 3);
rs = data(:, 4);
T  = data(:, 5);
Lp = data(:, 6);

eta     = visc(T);
lambdas = rs ./ rp;
parti   = (1 - lambdas) .^ 2;

D  = stokes_einstein(rs, eta, T);
Kc = kc_dech(lambdas);
Kd = kd_dech(lambdas);
Pe = peclet(rp, Kd, Kc, Jv, D, Lp);
Sm = sm_hindered(parti, Kc, Pe);
k  = kamicon8010(D, w, eta);

Sobs = sm2sobs(Sm, Jv, k);

