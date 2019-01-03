function [D, z, rg] = propRNA(str)
% PROPRNA  Returns properties of a RNA molecule identified by
% the string STR.
%   [D, Z, RG] = PROPRNA(STR) returns the diffusion coefficient
%   (D, [m2/s]), electric valence (Z) and radius of gyration
%   (RG, [m]) of a RNA molecule identified by the string STR.
%   The function accepts 3 different types of RNA molecules:
%     * RNA 23S
%     * RNA 16S
%     * RNA 5S
%   STR must be one of the following string:
%     * 'RNA23S'
%     * 'RNA16S'
%     * 'RNA5S'
%   The strings are not case sensitive. However white spaces
%   between characters should be avoided.
%   Radius of gyration and diffusion coefficients were obtained
%   in [1].
% 
%   References
%   
%   [1] J. C. Nunes et al, Modeling of Plasmid DNA/RNA separation
%       by ultrafiltration and application study (submitted).
%
%   see also propSAL

% Author: J. Nunes <josenunes34@gmail.com>
% Date: 2013

% String with function name
fn = 'propRNA';

if nargout~=3
    error(['Not enough output arguments\n'...
        'For further details type:\n',...
        '    help %s'], fn)
end

switch lower(str)
    case {'rna5s'}
        D  = 7.29e-11;
        z  = -120;
        rg = 1.505*3.39e-9;
        return
    case {'rna16s'}
        D  = 1.82e-11;
        z  = -1541;
        rg = 1.505*13.6e-9;
        return
    case {'rna23s'}
        D  = 1.39e-11;
        z  = -2904;
        rg = 1.505*17.8e-9;
        return
    otherwise
        error(['%s is not a valid RNA id.\n' ...
                'Use RNA23S or RNA16S or RNA5S.\n' ...
                'For further details type:\n',...
                '    help %s'], str, fn)
end