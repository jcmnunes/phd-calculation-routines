function parti = fiflex(lambdag)
% FIFLEX  evaluates the partition coefficient, for long 
% and flexible macromolecules, under dynamic conditions. 
%   PARTI = FIFLEX(LAMBDAG) evaluates the dynamic 
%   partition coefficient as a function of LAMBDAG (the
%   ratio of the radius of gyration to the membrane 
%   pore radius). The equation proposed by Morao et al [1]
%   is used. This equation is valid for both CSC and FJC
%   models. 
%
%   References
%
%   [1] Morao et al, Ultrafiltration of supercoiled plasmid
%       DNA: Modelling and application, Journal of Membrane 
%       Science 378 (2011) 280-289.

x  = log(lambdag) + 0.896;
a0 = -11.6;
a1 = 11.53;
a2 = 3.955;
a3 = -5.52;
a4 = -0.613;

parti = exp(a0 + a1 / (2 * a3) * (2 * a4 * ...
	log(exp((x + a3 / 2) / a4) + exp(a2 / a4)) - ...
	2 * a4 * log(exp((a2 + a3 / 2) / a4) + ...
	exp(x / a4)) + a3));
