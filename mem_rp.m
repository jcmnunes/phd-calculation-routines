function [rp_opt, graph] = mem_rp (data, rp1, rp2, RelTol)
% MEM_RP  Membrane characterization (pore radius) by
% the Hindered transport model [1].
%   RP_OPT = MEM_RP(DATA) evaluates the optimal RP value
%   that minimizes the Sy error function [2]. 
%   
%   There are many ways to use this function. Below, the
%   recommended procedure is described:
%   
%   1 - Import data to the workspace. DATA must be a 
%       cell array, and must be formatted as follows (without
%       column headers):
%       
%       Jv   w    T    Lp   Sobs  Solute 
%       --------------------------------
%       num  num  num  num  num   str  
%       num  num  num  num  num   str
%       ...
%
%       Jv     - filtration flux [m/s]
%       w      - stirring speed [rad/s]
%       T      - temperature [K]
%       Lp     - hydraulic permeability [m]
%       Sobs   - observed sieving coefficient
%       Solute - specifies the solute. Solute must be a linear
%                PEG or a linear dextran. Append to the solute
%                name the molecular weight (in Da).
%                For instance, if the solute is a PEG and
%                has a molecular weight equal to 20kDa,
%                Solute = 'peg 20000'.
%                If solute is a dextran, with Mw = 70 kDa,
%                Solute = 'dextran 70000'. For dextrans, it 
%                is also possible to use the 'T',
%                nomenclature. For instance: 
%                Solute = 'dex 70000' and Solute = 'dex T70'
%                both represent a linear dextran with 
%                Mw = 70 kDa. 
% 
%   2 - Call the function for the first time with only one
%       input argument (DATA), and make a plot of the Sy
%       function vs rp by providing a second output argument:
%       >> [rp_opt, graph] = mem_rp(data);
%
%   3 - Inspect the plot. See if the optimal rp value is
%       realistic and the function Sy is well behaved.
%       If not, it is possible to specify a new search 
%       interval. If the new interval is rp in [rp1 rp2] then:
%       >> [rp_opt, graph] = mem_rp(data, rp1, rp2);
%       will make a plot of Sy in the new interval.
%
%   4 - If you think the solution is acceptable, it is possible
%       to increase the precision, giving to the function a 
%       4th input argument (RelTol):
%       RelTol = 1e-5;
%       >> [rp_opt, graph] = mem_rp(data, rp1, rp2, RelTol);
%
%   If the above procedure doesn't work make sure DATA is 
%   formatted like the above template, and all units are 
%   SI units. Also make sure you entered the values correctly,
%   row by row. It is advisable to use more than one solute in 
%   membrane characterization. Also make sure all the values
%   in DATA were obtained with the same membrane, and the 
%   obtained observed sieving coefficients are neither
%   too close to one, or too close to zero.
%
%   References
%
%       [1] W. M. Deen, AIChE Journal 33 (1987) 1409-1425.
%       [2] W. R. Bowen et al, Journal of Membrane Science
%           126 (1997) 91-106. 
%
%   see also sm2sobs, kd_dech, kc_dech, peclet

% Author: J. Nunes <josenunes34@gmail.com>
% Date: 2013

% Initiate string with function name
fn = 'mem_rp';

% See if first input variable is of type cell
if ~iscell(data)
    error(['The first input variable must be of type cell\n', ...
        'For further details type:\n', ...
        '    help %s'], fn);
end

% See if data has 6 columns
if length(data(1, :)) ~= 6
    error(['The first input variable must have 6 columns.\n', ...
        'For further details type:\n', ...
        '    help %s'], fn);
end      

% Check if flux values are in SI units
x = [data{:, 1}];
if any(x > 1)
    error(['Flux values higher than 1 m/s detected.\n', ...
        'Make sure flux units are SI units (m/s).\n', ...
        'For further details type:\n', ...
        '    help %s'], fn)
end

% Temperature values in SI units?
x = [data{:, 3}];
if any(x < 200)
    error(['Temperatures lower than 200K detected.\n', ...
        'Make sure temperature units are SI units [K].\n', ...
        'For further details type:\n', ...
        '    help %s'], fn)
end

% Permeability values in SI units?
x = [data{:, 4}];
if any(x > 1)
    error(['Permeability values higher than 1 m detected.\n', ...
        'Make sure permeability units are SI units [m].\n', ...
        'For further details type:\n', ...
        '    help %s'], fn)
end

% Sobs higher then 1??
x = [data{:, 5}];
if any(x > 1)
    error(['Observed permeation values higher than 1 detected.\n',...
        'For further details type:\n', ...
        '    help %s'], fn)
end

% Sobs negative??
x = [data{:, 5}];
if any(x < 0)
    error(['Negative observed permeation values detected.\n',...
        'For further details type:\n', ...
        '    help %s'], fn)
end

% Check if rp2 > rp1
if nargin > 2
    if rp1 > rp2
       error(['rp1 > rp2 ...\n',...
        'For further details type:\n', ...
        '    help %s'], fn) 
    end
end

% Check rp1 and rp2
if nargin > 1
    if rp1 > 1e-3
        error(['rp1 higher than 1 mm... \n',...
        'For further details type:\n', ...
        '    help %s'], fn)
    end
    if rp1 < 0
        error(['rp1 negative... \n',...
        'For further details type:\n', ...
        '    help %s'], fn)
    end
end
if nargin > 2
    if rp2 > 1e-3
        error(['rp2 higher than 1 mm... \n',...
        'For further details type:\n', ...
        '    help %s'], fn)
    end
    if rp2 < 0
        error(['rp2 negative... \n',...
        'For further details type:\n', ...
        '    help %s'], fn)
    end
end

% Check RelTol 
if nargin > 3
    if RelTol > 1e-3
       error(['Use values lower or equal than 1e-3 for RelTol.\n',...
        'For further details type:\n', ...
        '    help %s'], fn) 
    end
    if RelTol < 1e-7
        error(['Use RelTol values in the interval [1e-3, 1e-7].\n',...
        'For further details type:\n', ...
        '    help %s'], fn) 
    end
end

tic

% Define the size of rp_range. Increase for better accuracy. 
ssize = 250;

% Allocate variables
len_data      = length(data(:,1));
deltaSobs     = ones(len_data, 1);
Sy            = ones(ssize, 1);
sumDeltaSobs  = ones(ssize, 1);

% Evaluate rs of solute from Solute string
rs = ones(len_data, 1);
for ii = 1:len_data
    str_SOLUTE = data{ii, 6};
    % Lower case
    str_solute = lower(str_SOLUTE);
    % Remove spaces
    str_solute(isspace(str_solute))=[];
    if str_solute(1) == 'd'
        pat = '\w\d+';
        solute_id = regexp(str_solute, pat, 'match');
        str_test = solute_id{1};
        if str_test(1) == 't'
            pat = '\d+';
            num_solute = regexp(solute_id, pat, 'match');
            Mw = str2double(num_solute{1}) * 1e3;
            rs(ii) = 0.0282 * Mw^0.47752 * 1e-9;
        else
            pat = '\d+';
            num_solute = regexp(solute_id, pat, 'match');
            Mw = str2double(num_solute{1});
            rs(ii) = 0.0282 * Mw^0.47752 * 1e-9;
        end
    elseif str_solute(1) == 'p'
        pat = '\d+';
        num_solute = regexp(str_solute, pat, 'match');
        Mw = str2double(num_solute{1});
        rs(ii) = exp((-23.91 + 0.4648 * log(Mw)));
    else
        error(['%s is not a valid solute id.\n'...
            'For further details type:\n'...
            '    help %s'], str_SOLUTE, fn)
    end
end

% Define the default values of rp1 and rp2
Sobs_exp = [data{:, 5}]';
c = mean(Sobs_exp);
if nargin < 2, rp1 = max(rs) - 0.2 * max(rs); end
if nargin < 3, rp2 = max(rs) + c * 1.5 * max(rs); end

% Range of rp values to find the minimum of Sy 
rp_range = linspace(rp1, rp2, ssize)';

% Execute 'for' cycle for every rp in rp_range
for ii = 1:ssize
    rp = rp_range(ii);
    for jj = 1:len_data
        
        Jv = data{jj, 1};
        w  = data{jj, 2};
        T  = data{jj, 3};
        Lp = data{jj, 4};
        
        eta     = visc(T);
        lambdas = rs(jj) / rp;
        D       = stokes_einstein(rs(jj), eta, T);
        Kc      = kc_dech(lambdas);
        Kd      = kd_dech(lambdas);
        parti   = (1 - lambdas)^2;
        Pe      = peclet(rp, Kd, Kc, Jv, D, Lp);
        k       = kamicon8010(D, w, eta);
        Sm      = sm_hindered(parti, Kc, Pe);
        
        Sobs_calc     = sm2sobs(Sm, Jv, k);
        deltaSobs(jj) = (Sobs_calc - Sobs_exp(jj))^2;
    
    end
    sumDeltaSobs(ii) = sum(deltaSobs);
    % Evaluate Sy function for rp_range(ii)
    Sy(ii) = sqrt(sumDeltaSobs(ii) / (len_data - 1));
end

% If requested higher precision fit a spline and find minimum 
if nargin > 3
    xx = linspace(rp1, rp2, 1 / RelTol);
    S = spline(rp_range, Sy, xx);
    [~, ind] = min(S);
    rp_opt = xx(ind);
else
    [~, index] = min(Sy);
    rp_opt = rp_range(index);
end

% plot section
str = sprintf('Optimal rp value is %s m', num2str(rp_opt));
if nargout > 1
    if nargin < 4
        figure('name', str)
        plot(rp_range, Sy, 'k-', rp_opt, min(Sy),'bo')
        xlim([rp1 rp2])
        xlabel('r_p [m]');
        ylabel('S_y');
        legend('S_y','Optimal r_p')
        graph = true;
    else
        figure('name', str)
        plot(xx, S, 'k-', rp_opt, min(S),'bo')
        xlim([rp1 rp2])
        xlabel('r_p [m]');
        ylabel('S_y');
        legend('S_y','Optimal r_p')
        graph = true;
    end
end
toc
disp(str)

