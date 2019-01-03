function [D1, z1, D2, z2] = propSAL(str)
% PROPSAL  Returns properties of a dissolved salt identified
% by the string STR.
%   [D1, Z1, D2, Z2] = PROPSAL(STR) returns the diffusion 
%   coefficient and electric valence of the salt anion
%   (D1 [m2/s] and Z1), and the diffusion coefficient
%   and electric valence of the salt cation (D2 [m2/s] and
%   Z2).
%   The function accepts 3 different salts:
%     * NaCl
%     * CH3COOK
%     * CaCl2
%   STR must be one of the following string:
%     * 'NaCl'
%     * 'CH3COOK'
%     * 'CaCl2'
%   The strings are not case sensitive. However white spaces
%   between characters should be avoided.
%   Diffusion coefficients were obtained in [1].
%
%   References
%   
%   [1] W. M. Haynes, CRC - Handbook of Chemistry and Physics. 
%
%   see also propRNA    

% Author: J. Nunes <josenunes34@gmail.com>
% Date: 2013


% String with function name
fn = 'propSAL';

if nargout~=4
    eerror(['Not enough output arguments\n'...
        'For further details type:\n',...
        '    help %s'], fn)
end

switch lower(str)
    case {'nacl'}
        D1=2.032e-9;D2=1.334e-9;
        z1=-1;z2=1;
        return
    case {'ch3cook'}
        D1=1.957e-9;D2=1.089e-9;
        z1=-1;z2=1;
        return
    case {'cacl2'}
        D1=2.032e-9;D2=0.792*2e-9;
        z1=-1;z2=2;
        return
    otherwise
        error(['%s is not a valid Salt id.\n' ...
                'Use NaCl or CH3COOK or CaCl2.\n' ...
                'For further details type:\n',...
                '    help %s'], str, fn)
end
end 