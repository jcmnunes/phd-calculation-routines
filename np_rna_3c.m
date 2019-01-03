function [Sobs, pro, logF, graph] = np_rna_3c (data, RelTol)
% NP_RNA_3C  Integration of the Nernst-Planck equation for an
% ionic system with 3 components (RNA and 2 salt ions) during 
% filtration in Amicon 8010 stirred cell.
%   SOBS = NP_RNA_3C(DATA) tries to integrate the system of
%   first-order differential equations defined by the Nernst-
%   -Planck equation applied to membrane filtration [1]. DATA 
%   must be a cell array containing both doubles and strings.
%   The system is defined for 3 components:
%     1. RNA (23S or 16S or 5S)
%     2. co-ion
%     3. counter-ion
%   The cell array DATA must be formatted as follows
%   (without column headers):
%     
%     Jv   I    w    rp   Cb1  T    RNA  Salt
%     ---------------------------------------
%     num  num  num  num  num  num  str  str
%     num  num  num  num  num  num  str  str
%     ...
% 
%     Jv   - filtration flux [m/s]
%     I    - ionic strength [mol/m3]
%     w    - stirring speed [rad/s]
%     rp   - membrane pore radius [m]
%     Cb1  - RNA concentration in feed solution [mol/m3]
%     T    - temperature [K]
%     RNA  - specifies the RNA species. Must be one of
%            the following strings:
%              * RNA23S
%              * RNA16S
%              * RNA5S
%     Salt - specifies the salt. Must be one of the 
%            following strings:
%              * NaCl
%              * CH3COOK
%              * CaCl2
%     These strings are not case sensitive. However, spaces between
%     characters should be avoided.  
%          
%   SI units must be used.
%   
%   SOBS is the calculated observed permeation for the
%   RNA species. 
%     
%   Example:
%   % Initiate the cell array DATA (note the use of {} delimiters)
%   data = {1e-6, 150, 760*2*pi/60, 10e-9, 1e-5, 298, 'RNA23S', ...
%           'NaCl'};
%   % call function
%   Sobs = np_rna_3c(data)
%   
%   Sobs =
%       
%       0.2052
%         
%   Sometimes, however, this method may not be the most
%   appropriate, especially for those cases when the 
%   number of rows in cell array DATA is high. If that is
%   the case, it is better to use a spreadsheet to organize
%   the input data, formatting the tables according to the
%   above template. The data can then be easily imported
%   in a cell array format (for instance, using the XLSREAD
%   function). For MATLAB users the Spreadsheet Link(TM) EX
%   add-in may be a better alternative. For GNU Octave users
%   the IO package is required to run the XLSREAD function
%   (see <http://octave.sourceforge.net/>).
%     
%   [SOBS, PRO] = NP_RNA_3C(DATA) PRO is a structure array
%   containing the result of the integration. PRO contains
%   5 fields:
%     * dist      - the spatial coordinate used in the 
%                   integration
%     * rna       - the RNA concentration profile
%     * anion     - the salt anion concentration profile
%     * cation    - the salt cation concentration profile
%     * potential - the electric potential profile
%   
%   Example: Access the salt anion concentration profile and
%   plot it along the spatial coordinate, for the 5th row 
%   values in DATA.
%     
%     - Extract the salt anion concentration profile
%       from PRO structure: 
%         >> c = PRO(5).anion;
%     - Extract the spatial coordinate
%         >> y = PRO(5).dist;
%     - Plot the result:
%         >> plot(y, c)
%  
%   [SOBS, PRO, LOGF] = NP_RNA_3C(DATA) LOGF is a cell array
%   with the following entries:
%          
%     logF =
% 
%     testODEflag      Cm1e  Cm1c  Cm2e  Cm2c
%     ----------------------------------------
%     OK (or failed)   num   num   num   num
%     ...
%           
%   testODEflag returns the string 'OK' if the solution 
%   obtained passes the ISCRES function test. This test
%   is based on the observation that some solutions are 
%   bad, because the integration method fails. This is
%   a first indicator of the quality of the solution.
%   * Cm1e is the last estimated value for the RNA 
%        concentration at the membrane surface.
%   * Cm1c is the last calculated value of the RNA
%        concentration near the membrane surface. 
%   * Cm2e is the last estimated value for the salt anion
%        concentration at the membrane surface.
%   * Cm2c is the last calculated value of the salt anion
%        concentration near the membrane surface.
%   If everything goes as expected, each of these pairs of
%   values (estimated and calculated) should be identical,
%   with a relative error defined by RelTol (see below).   
%   The variable LOGF thus offers a fast method to 
%   evaluate the quality of the solutions.
%     
%   [SOBS, PRO, LOGF, GRAPH] = NP_RNA_3C(DATA) if a 4th output
%   variable is specified, the function will make a plot
%   of the result of the integration, for every row in DATA.
%
%   As usual, it is possible to use the placeholder '~' to 
%   specify the output variables' positions.
%   
%     Example: Use the function just to plot the results of
%     the integration.
%     
%       >> [~, ~, ~, GRAPH] = NP_RNA_3C(DATA)
% 
%     
%   It is possible to change the default relative tolerance
%   (RelTol) of the function, defined as 1e-4. Just specify
%   the new RelTol as the second input argument. RelTol is
%   calculated by the formula:
%     
%     >> abs(Cme - Cmc) / min(Cme, Cmc)
% 
%   where Cme is the estimated value, Cmc is the calculated
%   value and ABS and MIN return the absolute and minimum 
%   values, respectively. 
%     
%     Example: Define a new RelTol
%     
%       >> RelTol = 1e-5;
%       >> Sobs   = np_rna_3c(data, RelTol)
% 
%   Sometimes, the function can't find a solution. This
%   is especially true for SOBS values close to 1. When 
%   this happens, the function issues a warning message
%   and points to the problematic row(s) within DATA.
%
%   For more information on the Nernst-Planck ODE system
%   see [1].
%
%   References
%    
%   [1] Morao et al, Journal of Membrane Science 378 
%       (2011) 280-289. 
%  
%   See also XLSREAD, ISCRES, NERNST3C.

% Author: J. Nunes <josenunes34@gmail.com>
% Date: 2013

% Initiate string with function name
fn = 'np_rna_3c';

% See if first input variable is of type cell
if ~iscell(data)
  error(['The first input variable must be of type cell\n', ...
         'For further details type:\n', ...
         '    help %s'], fn)
end

% See if data has 8 columns
if length(data(1, :)) ~= 8
  error(['The first input variable must have 8 columns.\n', ...
         'For further details type:\n', ...
         '    help %s'], fn)
end    

% Check if flux values are in SI units
x = [data{:, 1}];
if any(x > 1)
  error(['Flux values higher than 1 m/s detected.\n', ...
         'Make sure flux units are SI units (m/s).\n', ...
         'For further details type:\n', ...
         '    help %s'], fn)
end

% See if flux values are unusually high
if any(x > 1e-3)
  warning(['Flux values higher than 1 mm/s detected.\n', ...
           'Function may not find a solution.\n', ...
           'Make sure flux units are SI units (m/s).\n', ...
           'Use ctrl+c to stop the function if needed.\n', ...
           'For further details type:\n', ...
           '    help %s'], fn)
end

% Membrane pore radius values in SI units?
x = [data{:, 4}];
if any( x > 1e-3)
  error(['Membrane pore radius values higher than 1 mm detected.\n', ...
         'Make sure pore radius units are SI units.\n', ...
         'For further details type:\n', ...
         '    help %s'], fn)
end

% Cb1 values in SI units? (Cb1 - RNA concentration)
x = [data{:, 5}];
if any(x > 1)
  error(['Unusually large RNA concentration values.\n', ...
         'Make sure RNA concentration units are SI units.\n', ...
         'For further details type:\n', ...
         '    help %s'], fn)
end

% Temperature values in SI units?
x = [data{:, 6}];
if any(x < 200)
  error(['Temperatures lower than 200K detected.\n', ...
         'Make sure temperature units are SI units.\n', ...
         'For further details type:\n', ...
         '    help %s'], fn)
end

% Check RelTol 
if nargin > 1
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

% Define the default relative tolerance
if nargin < 2, RelTol = 1e-4; end

% Number of rows in data
len = length(data(:, 1));

% Allocate output variables
Sobs = zeros(len, 1);

if nargout >= 2
  pro(len).dist      = [];
  pro(len).rna       = [];
  pro(len).anion     = [];
  pro(len).cation    = [];
  pro(len).potential = [];
end

if nargout >= 3
  logF = {'testODEflag','Cm1 est','Cm1 calc','Cm2 est','Cm2 calc'};
end
% The 4th output variable does not need pre-allocation (plots)


% MAIN CYCLE
for jj = 1 : len
  counter    = jj;
  strcounter = num2str(counter);
  % SI units
  Jv      = data{jj, 1};
  I       = data{jj, 2};
  w       = data{jj, 3};
  rp      = data{jj, 4};
  Cb1     = data{jj, 5};
  T       = data{jj, 6};
  strRNA  = data{jj, 7};
  strSALT = data{jj, 8};
    
  [D1, z1, rg]     = propRNA(strRNA);
  [D2, z2, D3, z3] = propSAL(strSALT);
  F                = 96485.3365; % Faraday's constant
  R                = 8.3144621;  % Ideal gas constant
    
  % Calculate anion concentration from ionic strength
  switch lower(strSALT)
    case {'nacl'}
      Cb2    = I;
      anion  = 'Cl^-';
      cation = 'Na^+';
    case {'ch3cook'}
      Cb2    = I;
      anion  = 'CH_3COO^-';
      cation = 'K^+';
    case {'cacl2'}
      Cb2    = I / 2;
      anion  = 'Cl^-';
      cation = 'Ca^{2+}';
    otherwise
      error(['%s is not a valid salt id.\n' ...
             'Use NaCl or CH3COOK or CaCl2.\n' ...
             'For further details type:\n', ...
             '    help %s'], strSALT, fn)
  end
  Cb3  = - (z1 * Cb1 + z2 * Cb2) / z3;
  Psi0 = 0;
  % Initial conditions for ODE solver (rk4)
  x0 = [Cb1; Cb2; Cb3; Psi0];
    
  visco   = visc(T);
  lambdag = rg / rp;
  fi      = fiflex(lambdag);
  k       = kamicon8010(D1, w, visco);
  delta   = D1 / k;
    
  Cm1  = Cb1;
  Cm2  = Cb2;
  cond = 0;
    
  % MAIN 'while' CYCLE
  while(cond == 0)
    sign1 = 1; sign2 = 1;
    % Cm1 optimization
    while sign1 == sign2
      [Cm1, Cm2, condout] = testODE_3c(Cm1, Cm2, fi, ...
        delta, x0, z1, z2, z3, F, R, T, Jv, D1, D2, D3);
      if condout == 1, break, end
      Cp1 = fi * Cm1; Cp2 = Cm2;
      Cp3 = - (z1 * Cp1 + z2 * Cp2) / z3;
      [~, f] = rk4arg(@nernst3c, [0 delta], x0, 100, ...
        Cp1, Cp2, Cp3, z1, z2, z3, F, R, T, Jv, D1, ...
        D2, D3);
      dif1 = Cm1 - f(end, 1); if dif1 == 0, break, end
      Cm1e1 = Cm1;
      [mnt, expo] = mantexpnt(dif1);
      Cm1e2 = Cm1 - (sign(mnt) * ceil(abs(mnt))) * 10 ^ (expo);
      Cm1 = Cm1e2;
      [Cm1, Cm2, condout] = testODE_3c(Cm1, Cm2, fi, ...
        delta, x0, z1, z2, z3, F, R, T, Jv, D1, D2, D3);
      if condout == 1, break, end
      Cm1e2 = Cm1;
      Cp1 = fi * Cm1; Cp2 = Cm2;
      Cp3 = - (z1 * Cp1 + z2 * Cp2) / z3;
      [~, f] = rk4arg(@nernst3c, [0 delta], x0, 100, ...
        Cp1, Cp2, Cp3, z1, z2, z3, F, R, T, Jv, D1, ...
        D2, D3);
      dif2 = Cm1 - f(end, 1); if dif2 == 0, break ,end
      sign1 = sign(dif1); sign2 = sign(dif2);
      [mnt, expo] = mantexpnt(dif2);
      Cm1 = Cm1 - (sign(mnt) * ceil(abs(mnt))) * 10 ^ (expo);
    end
    if condout == 1, break, end
    Cm1 = Cm1e1 - (Cm1e2 - Cm1e1) / (dif2 - dif1) * dif1;
    sign11 = 1; sign22 = 1;
        
    % Cm2 optimization
    while sign11 == sign22
      Cp1 = fi * Cm1; Cp2 = Cm2;
      Cp3 = - (z1 * Cp1 + z2 * Cp2) / z3;
      [~, f] = rk4arg(@nernst3c, [0 delta], x0, 100,...
          Cp1, Cp2, Cp3, z1, z2, z3, F, R, T, Jv, D1, ...
          D2, D3);
      dif11 = Cm2 - f(end, 2); if dif11 == 0, break, end
      Cm2e1 = Cm2;
      [mnt, expo] = mantexpnt(dif11);
      Cm2e2 = Cm2 - (sign(mnt) * ceil(abs(mnt))) * 10 ^ (expo);
      Cm2 = Cm2e2; Cp2 = Cm2;
      Cp3 = - (z1 * Cp1 + z2 * Cp2) / z3;
      [y, f] = rk4arg(@nernst3c, [0 delta], x0, 100,...
          Cp1, Cp2, Cp3, z1, z2, z3, F, R, T, Jv, D1, ...
          D2, D3);
      dif22 = Cm2 - f(end, 2); if dif22 == 0, break, end
      % Variables for RelTol test
      d1 = abs(Cm1 - f(end, 1)); m1 = min(Cm1, f(end, 1));
      d2 = abs(Cm2 - f(end, 2)); m2 = min(Cm2, f(end, 2));
      % -----------------------------------
      sign11 = sign(dif11); sign22 = sign(dif22);
      [mnt, expo] = mantexpnt(dif22);
      Cm2 = Cm2 - (sign(mnt) * ceil(abs(mnt))) * 10 ^ (expo);
    end
    Cm2 = Cm2e1 - (Cm2e2 - Cm2e1) / (dif22 - dif11) * dif11;
        
    % RelTol test
    tol1 = d1 / m1; tol2 = d2 / m2;
    if tol1 < RelTol
      if tol2 < RelTol
        cond = 1;
      else
        cond = 0;
      end
    else
      cond = 0;
    end
  end % End of main 'while' cycle
    
  % Output variables for DATA(JJ, :)
    
  % SOBS (RNA observed sieving coefficient)
  if condout == 1
    Sobs(jj) = NaN;
    fprintf(['\nWarning:\nFunction did not find a\n'...
             'solution for row %s of input cell.\n'...
             'Please recheck the input values.\n'...
             'Also, Sobs may be too close to 1.\n']...
             ,strcounter)
  else
    Sobs(jj) = fi * Cm1 / Cb1;
  end
    
  % PRO structure
  if nargout >= 2
    if condout == 0
      pro(jj).dist      = y;
      pro(jj).rna       = f(:, 1);
      pro(jj).anion     = f(:, 2);
      pro(jj).cation    = f(:, 3);
      pro(jj).potential = f(:, 4);
    else
      pro(jj).dist      = 'No solution found';
      pro(jj).rna       = 'No solution found';
      pro(jj).anion     = 'No solution found';
      pro(jj).cation    = 'No solution found';
      pro(jj).potential = 'No solution found';
    end
  end
    
  % LOGF cell
  if nargout >= 3
    if condout == 1
      logF{jj+1, 1} = 'No solution found';
      logF{jj+1, 2} = 'No solution found';
      logF{jj+1, 3} = 'No solution found';
      logF{jj+1, 4} = 'No solution found';
      logF{jj+1, 5} = 'No solution found';
    else
      finaltest = iscres(f(:, 1));
      if finaltest == 1
        flag = 'OK!';
      else
        flag = 'Failed!';
      end    
        logF{jj+1, 1} = flag;
        logF{jj+1, 2} = Cm1;
        logF{jj+1, 3} = f(end, 1);
        logF{jj+1, 4} = Cm2e2;
        logF{jj+1, 5} = f(end, 2);
    end
  end
    
  % Plot polarization profiles
  % New figure for each row in data
  if nargout >= 4
    if condout == 1
      figstr = ['No solution found for values in row ', ...
                 num2str(counter)];
      figure('name',figstr)
      graph = true;
    else
      figstr = ['Polarization for values in row ', ...
                 num2str(counter),' of input data'];
      figure('name', figstr)
      subplot(2, 2, 1)
      plot(y, f(:, 1), 'k-')
      title(strRNA)
      xlim([0 delta])
      xlabel('y [m]')
      ylabel('C_1 [mol/m^3]')
      subplot(2, 2, 2)
      plot(y, f(:, 2), 'k-')
      title(anion)
      xlim([0 delta])
      xlabel('y [m]')
      ylabel('C_2 [mol/m^3]')
      subplot(2, 2, 3)
      plot(y, f(:, 3), 'k-')
      title(cation)
      xlim([0 delta])
      xlabel('y [m]')
      ylabel('C_3 [mol/m^3]')
      subplot(2, 2, 4)
      plot(y, f(:, 4), 'k-')
      title('potential')
      xlim([0 delta])
      xlabel('y [m]')
      ylabel('\Psi [V]')
      graph = true;
    end
  end
    
end % End of main 'for' cycle

strtoc = '\nElapsed time: %.6f seconds\n';
fprintf(strtoc, toc)