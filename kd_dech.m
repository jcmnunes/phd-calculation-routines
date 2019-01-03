function Kd = kd_dech(l)
% KD_DECH  Implementation of Dechadilok and Deen's
% equation to compute the hindrance diffusion
% coefficient [1].
% 	KD = KD_DECH(L) computes de hindrance diffusion
%   coefficient. L is the ratio of solute Stokes radius 
%   to the pore radius. 
%   For more information see [1].
%
%   References
%
%   [1] P. Dechadilok and W. M . Deen, Industrial and 
%       Engineering Chemistry Research, 45 (2006) 
%       6953-6959. 
%
%   see also kc_dech

fi = (1 - l).^2;

a = -1.56034;
b = 0.528155;
c = 1.91521;
d = -2.81903;
e = 0.270788;
f = 1.10115;
g = -0.435933;

Kd = (1 + 9 .* l .* log(l) / 8 + a * l + b * l.^2 + c * ...
	l.^3 + d * l.^4 + e * l.^5 + f * l.^6 ...
  + g * l.^7) ./ fi;