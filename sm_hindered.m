function Sm = sm_hindered(parti, Kc, Pe)
% SM_HINDERED  Intrinsic permeation by the
% hindered transport model.
%   SM = SM_HINDERED(PARTI, KC, PE) evaluates the
%   intrinsic permeation coefficient (SM) by the
%   hindered transport model. 
%   Input variables are (adimensional):
%     * PARTI - Partition coefficient
%     * KC    - Hindered convection coefficient
%     * PE    - Peclet number

Sm = parti .* Kc ./ (1 - (1 - parti .* Kc) .* exp(-Pe));