function rg=rgcsc(nk,L)
% RGCSC  Compute the radius of gyration of a CSC.
%   RG = RGCSC(NK, L) returns the radius of 
%   gyration (RG, [m]) of a CSC, as a function
%   of the number of segments (NK) and the chain
%   length (L, [m]). The correlation proposed by
%   Morao et al is used [1]. 
% 
%   References
% 
%   [1] Morao et al, Ultrafiltration of supercoiled plasmid
%       DNA: Modelling and application, Journal of Membrane 
%       Science 378 (2011) 280-289.

a = 0.221;
b = -0.352;

rg = a * nk ^ b * L;