function Kc = kc_dech(l)
% KC_DECH  Implementation of the Dechadilok and Deen's
% equation to compute the hindrance convection
% coefficient [1].
% 	KC = kC_DECH(L) computes de hindrance convection
%   coefficient. L is the ratio of solute radius 
%   to the pore radius. 
%   For more information see [1].
%
%   References
%
%   [1] P. Dechadilok and W. M . Deen, Industrial and 
%       Engineering Chemistry Research, 45 (2006) 
%       6953-6959. 
%
%   see also kd_dech

a = 3.867;
b = -1.907;
c = -0.834;
d = 1.867;
e = -0.741;

Kc = (1 + a * l + b * l.^2 + c * l.^3)./ ... 
     (1 + d * l + e * l.^2);

