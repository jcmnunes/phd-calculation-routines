function eta = visc(T)
% VISC  Viscosity of water
%   ETA = VISC(T) returns pure water dynamic viscosity (ETA, [Pa.s])
%   as a function of temperature (T, [K]).
%   The correlation of Morao is used [1].
%   The correlation was obtained from the values of Lide [2].
% 
%   References
% 
%   [1] A. M. Morao, Remocao de poluentes persistentes por
%       tecnologias de membranas e processos eletroquimicos,
%       Tese de Doutoramento, Universidade da Beira Interior,
%       2006.
% 
%   [2] Lide, D.R. (Ed.), Handbook of Chemistry and Physics, 
%       72nd Edition, CRC Press, Boca Raton.
    
eta = (0.398 + 1.39 .* exp(-(T - 273.15) / 23.8)) * 0.001;