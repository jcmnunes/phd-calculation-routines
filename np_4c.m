function [Sobs1, Sobs2, pro, logF, graph] = np_4c (data, RelTol)
% NP_4C  Integration of the Nernst-Planck equation for an
% ionic system with 4 components (pDNA(sc), RNA and 2 
% salt ions) during filtration in Amicon 8010 stirred cell.
%   [SOBS1, SOBS2] = NP_4C(DATA) tries to integrate the 
%   system of first-order differential equations defined
%   by the Nernst-Planck equation applied to membrane
%   filtration [1]. DATA must be a cell array containing
%   both doubles and strings. The system is defined for
%   4 components:
%     1. pDNA (sc)
%     2. RNA (23S or 16S or 5S)
%     3. co-ion
%     4. counter-ion
%   The cell array DATA must be formatted as follows
%   (without column headers):
%     
%     Jv   I    w    rp   Cb1  Cb2  nbp  T    RNA  Salt
%     -------------------------------------------------
%     num  num  num  num  num  num  num  num  str  str
%     num  num  num  num  num  num  num  num  str  str
%     ...
%
%     Jv   - filtration flux [m/s]
%     I    - ionic strength [mol/m3]
%     w    - stirring speed [rad/s]
%     rp   - membrane pore radius [m]
%     Cb1  - pDNA concentration in feed solution [mol/m3]
%     Cb2  - RNA concentration in feed solution [mol/m3]
%     nbp  - number of base pairs (pDNA)
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
%   These strings are not case sensitive. However, spaces between
%   characters should be avoided.
%          
%   SI units must be used.
%     
%   SOBS1 is the calculated observed permeation for the pDNA.
%   Sobs2 is the calculated observed permeation for RNA.
%
%   Example:
%   % Initiate the cell array DATA (note the use of {} delimiters)
%   data = {2e-6, 316, 2*pi*760/60, 15e-9, 4.47e-6, 2.26e-3, ...
%           6050, 298, 'rna5s', 'ch3cook'};
%   % call function
%   [Sobs1, Sobs2] = np_4c(data)
%     
%     Sobs1 =
%  
%         0.0911
% 
%     Sobs2 =
% 
%         0.8245
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
%   [SOBS1, Sobs2, PRO] = NP_4C(DATA) PRO is a structure 
%   array containing the result of the integration. 
%   PRO contains 5 fields:
%     * dist      - the spatial coordinate used in the 
%                   integration
%     * pdna      - the pDNA concentration profile       
%     * rna       - the RNA concentration profile
%     * anion     - the salt anion concentration profile
%     * cation    - the salt cation concentration profile
%     * potential - the electric potential profile
%     
%   Example: Access the salt cation concentration profile and
%   plot it along the spatial coordinate, for the 3rd row val-
%   ues in DATA.
%     
%     - First extract the sal cation concentration profile
%       from PRO structure: 
%         >> c = PRO(3).cation;
%     - Extract the spatial coordinate
%         >> y = PRO(3).dist;
%     - Plot the result:
%         >> plot(y, c)
%
%   [SOBS1, SOBS2, PRO, LOGF] = NP_4C(DATA) LOGF is a cell
%   array with the following entries:
%          
%     logF = 
% 
%       testODEflag      Cm1e  Cm1c  Cm2e  Cm2c  Cm3e  Cm3c 
%       ----------------------------------------------------
%       OK (or failed)   num   num   num   num   num   num
%       ...
%           
%   testODEflag returns the string 'OK' if the solution 
%   obtained passes the ISCRES function test. This test
%   is based on the observation that some solutions are 
%   bad, because the integration method fails. This is
%   a first indicator of the quality of the solution.
%   * Cm1e is the last estimated value for the pDNA 
%        concentration at the membrane surface.
%   * Cm1c is the last calculated value of the pDNA 
%        concentration near the membrane surface. 
%   * Cm2e is the last estimated value for the RNA 
%        concentration at the membrane surface.
%   * Cm2c is the last calculated value of the RNA 
%        concentration near the membrane surface.
%   * Cm3e is the last estimated value for the salt anion
%        concentration at the membrane surface.
%   * Cm3c is the last calculated value of the salt anion
%        concentration near the membrane surface.
%   If everything goes as expected, each of these pairs of
%   values (estimated and calculated) should be identical,
%   with a relative error defined by RelTol (see below).   
%   The variable LOGF thus offers a fast method to 
%   evaluate the quality of the solutions.
%
%   [SOBS1, Sobs2, PRO, LOGF, GRAPH] = NP_4C(DATA) if a 
%   5th output variable is specified, the function will
%   make a plot of the result of the integration, for every
%   row in DATA.
% 
%   As usual, it is possible to use the placeholder '~' to 
%   specify the output variables' positions.
%     
%     Example: Use the function just to plot the results of
%     the integration.
%     
%       >> [~, ~, ~, ~, GRAPH] = NP_4C(DATA)
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
%       >> [Sobs1, Sobs2] = np_4c(data, RelTol)
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
%   See also XLSREAD, ISCRES, NERNST4C.

% Author: J. Nunes <josenunes34@gmail.com>
% Date: 2013

% Initiate string with function name
fn = 'np_4c';

% See if first input variable is of type cell
if ~iscell(data)
  error(['The first input variable must be of type cell\n', ...
         'For further details type:\n', ...
         '    help %s'], fn)
end     

% See if data has 10 columns
if length(data(1, :)) ~= 10
  error(['The first input variable must have 10 columns.\n', ...
         'For further details type:\n', ...
         '    help %s'], fn)
end 

% Check if flux is in SI units
x=[data{:, 1}];
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

% See if membrane pore radius values are in SI units
x = [data{:, 4}];
if any( x > 1e-3)
  error(['Membrane pore radius values higher than 1 mm detected.\n', ...
         'Make sure pore radius units are SI units.\n', ...
         'For further details type:\n', ...
         '    help %s'], fn)
end

% Cb1 values in SI units? (Cb1 - pDNA concentration)
x = [data{:, 5}];
if any(x > 1)
  error(['Unusually large pDNA concentration values.\n', ...
         'Make sure pDNA concentration units are SI units.\n', ...
         'For further details type:\n', ...
         '    help %s'], fn)
end

% Cb2 values in SI units? (Cb2 - RNA concentration)
x = [data{:, 6}];
if any(x > 1)
  error(['Unusually large RNA concentration values.\n',...
         'Make sure RNA concentration units are SI units.\n', ...
         'For further details type:\n\n',...
         ' help %s'],fn)
end

% Temperature values in SI units?
x = [data{:, 8}];
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

% Number of rows in DATA
len = length(data(:, 1));

% Allocate output variables
Sobs1 = zeros(len, 1);

if nargout >= 2
  Sobs2 = zeros(len, 1);
end

if nargout >= 3
  pro(len).dist      = [];
  pro(len).pdna      = [];
  pro(len).rna       = [];
  pro(len).anion     = [];
  pro(len).cation    = [];
  pro(len).potential = [];
end

if nargout >= 4
    logF = {'testODEflag','Cm1 est','Cm1 calc',...
        'Cm2 est','Cm2 calc','Cm3 est','Cm3 calc'};
end
% The 5th output variable does not need allocation (plots)


% MAIN CYCLE
for jj = 1 : len 
  counter = jj;
  strcounter = num2str(counter);
  % SI units
  Jv      = data{jj, 1};
  I       = data{jj, 2};
  w       = data{jj, 3};
  rp      = data{jj, 4};
  Cb1     = data{jj, 5};
  Cb2     = data{jj, 6};
  npb     = data{jj, 7};
  T       = data{jj, 8};
  strRNA  = data{jj, 9};
  strSALT = data{jj, 10};
    
  [D2, z2, rg2]    = propRNA(strRNA);
  lambdag2         = rg2 / rp;
  [D3, z3, D4, z4] = propSAL(strSALT);
  z1               = -2 * npb;
    
  F = 96485.3365; % Faraday's constant
  R = 8.3144621;  % Ideal gas constant
    
  % Calculate anion concentration from ionic strength
  switch lower(strSALT)
    case {'nacl'}
      Cb3    = I;
      cation = 'Na^+';
      anion  = 'Cl^-';
    case {'ch3cook'}
      Cb3    = I;
      cation = 'K^+';
      anion  = 'CH_3COO^-';
    case {'cacl2'}
      Cb3    = I/2;
      cation = 'Ca^{2+}';
      anion  = 'Cl^-';
    otherwise
      error(['%s is not a valid salt id.\n' ...
             'Use NaCl or CH3COOK or CaCl2.\n' ...
             'For further details type:\n',...
             '    help %s'], strSAL, fn)
  end
  
  Cb4  = - (z1 * Cb1 + z2 * Cb2 + z3 * Cb3) / z4;
  Psi0 = 0;
  % Initial conditions for ODE solver (rk4) 
  x0 = [Cb1; Cb2; Cb3; Cb4; Psi0];
    
  visco    = visc(T);
  lK       = manning(I, z4);
  lk       = 0.196 * lK;
  L        = npb * 0.34e-9;
  nk       = L / lk;
  rg1      = rgcsc(nk, L);
  lambdag1 = rg1 / rp;
  fi1      = fiflex(lambdag1);
  fi2      = fiflex(lambdag2);
  D1       = dpdnapra(T, visco, npb);
  k        = kamicon8010(D1, w, visco);
  delta    = D1 / k;
    
  Cm1  = Cb1;
  Cm2  = Cb2;
  Cm3  = Cb3;
  cond = 0;
    
  % MAIN 'while' CYCLE
  while(cond==0)
    sign1=1;sign2=1;    
    % Cm1 optimization
    while sign1==sign2
      [Cm1, Cm2, Cm3, condout] = testODE_4c(Cm1, Cm2, Cm3, fi1,...
        fi2, delta, x0, z1, z2, z3, z4, F, R, T, Jv, D1, D2, D3, D4);
      if condout==1, break, end
      Cp1 = fi1 * Cm1; Cp2 = fi2 * Cm2; Cp3 = Cm3;
      Cp4 = - (z1 * Cp1 + z2 * Cp2 + z3 * Cp3) / z4;
      [~, f] = rk4arg(@nernst4c, [0 delta], x0, 100, Cp1, Cp2, ...
        Cp3, Cp4, z1, z2, z3, z4, F, R, T, Jv, D1, D2, D3, D4);
      dif1 = Cm1 - f(end, 1); if dif1 == 0, break, end
      Cm1e1 = Cm1;
      [mnt, expo] = mantexpnt(dif1);
      Cm1e2 = Cm1 - (sign(mnt) * ceil(abs(mnt))) * 10 ^ (expo);
      Cm1 = Cm1e2;
      [Cm1, Cm2, Cm3, condout] = testODE_4c(Cm1, Cm2, Cm3, fi1,...
        fi2, delta, x0, z1, z2, z3, z4, F, R, T, Jv, D1, D2, D3, D4);
      if condout==1, break, end
      Cm1e2 = Cm1;
      Cp1 = fi1 * Cm1; Cp2 = fi2 * Cm2; Cp3 = Cm3;
      Cp4 = - (z1 * Cp1 + z2 * Cp2 + z3 * Cp3) / z4;
      [~, f] = rk4arg(@nernst4c, [0 delta], x0, 100, Cp1, Cp2, ...
        Cp3, Cp4, z1, z2, z3, z4, F, R, T, Jv, D1, D2, D3, D4);
      dif2 = Cm1 - f(end, 1); if dif2 == 0, break, end
      sign1 = sign(dif1); sign2 = sign(dif2);
      [mnt, expo] = mantexpnt(dif2);
      Cm1 = Cm1 - (sign(mnt) * ceil(abs(mnt))) * 10 ^ (expo);
    end
    if condout == 1, break, end
    Cm1 = Cm1e1 - (Cm1e2 - Cm1e1) / (dif2 - dif1) * dif1;
    sign11 = 1; sign22 = 1;
        
    % Cm2 optimization
    while sign11 == sign22
      Cp1 = fi1 * Cm1; Cp2 = fi2 * Cm2; Cp3 = Cm3;
      Cp4 = - (z1 * Cp1 + z2 * Cp2 + z3 * Cp3) / z4;
      [~, f] = rk4arg(@nernst4c, [0 delta], x0, 100, Cp1, Cp2, ...
        Cp3, Cp4, z1, z2, z3, z4, F, R, T, Jv, D1, D2, D3, D4);
      dif11 = Cm2 - f(end, 2); if dif11 == 0, break, end
      Cm2e1 = Cm2;
      [mnt, expo] = mantexpnt(dif11);
      Cm2e2 = Cm2 - (sign(mnt) * ceil(abs(mnt))) * 10 ^ (expo);
      Cm2 = Cm2e2; Cp2 = fi2 * Cm2;
      Cp4 = - (z1 * Cp1 + z2 * Cp2 + z3 * Cp3) / z4;
      [~, f] = rk4arg(@nernst4c, [0 delta], x0, 100, Cp1, Cp2, ...
        Cp3, Cp4, z1, z2, z3, z4, F, R, T, Jv, D1, D2, D3, D4);
      dif22 = Cm2 - f(end, 2); if dif22 == 0, break, end
      sign11 = sign(dif11); sign22 = sign(dif22);
      [mnt, expo] = mantexpnt(dif22);
      Cm2 = Cm2 - (sign(mnt) * ceil(abs(mnt))) * 10 ^ (expo);
    end
    Cm2=Cm2e1-(Cm2e2-Cm2e1)/(dif22-dif11)*dif11;
    sign111=1;sign222=1;
        
    % Cm3 optimization
    while sign111 == sign222
      Cp1 = fi1 * Cm1; Cp2 = fi2 * Cm2; Cp3 = Cm3;
      Cp4 = - (z1 * Cp1 + z2 * Cp2 + z3 * Cp3) / z4;
      [~, f] = rk4arg(@nernst4c, [0 delta], x0, 100, Cp1, Cp2, ...
        Cp3, Cp4, z1, z2, z3, z4, F, R, T, Jv, D1, D2, D3, D4);
      dif111 = Cm3 - f(end, 3); if dif111 == 0, break, end
      Cm3e1 = Cm3;
      [mnt, expo] = mantexpnt(dif111);
      Cm3e2 = Cm3 - (sign(mnt) * ceil(abs(mnt))) * 10 ^ (expo);
      Cm3 = Cm3e2; Cp3 = Cm3;
      Cp4 = - (z1 * Cp1 + z2 * Cp2 + z3 * Cp3) / z4;
      [y, f] = rk4arg(@nernst4c, [0 delta], x0, 100, Cp1, Cp2, ...
        Cp3, Cp4, z1, z2, z3, z4, F, R, T, Jv, D1, D2, D3, D4);
      % Variables for RelTol test
      d1 = abs(Cm1 - f(end, 1)); m1 = min(Cm1, f(end, 1));
      d2 = abs(Cm2 - f(end, 2)); m2 = min(Cm2, f(end, 2));
      d3 = abs(Cm3 - f(end, 3)); m3 = min(Cm3, f(end, 3));
      % -----------------------------------------
      dif222 = Cm3 - f(end, 3); if dif222 == 0, break, end
      sign111 = sign(dif111); sign222 = sign(dif222);
      [mnt, expo] = mantexpnt(dif222);
      Cm3 = Cm3 - (sign(mnt) * ceil(abs(mnt))) * 10 ^ (expo);
    end
    Cm3 = Cm3e1 - (Cm3e2 - Cm3e1) / (dif222 - dif111) * dif111;
        
    % RelTol test
    tol1 = d1 / m1; tol2 = d2 / m2; tol3 = d3 / m3;
    if tol1 < RelTol
      if tol2 < RelTol
        if tol3 < RelTol
          cond = 1;
        else
          cond=0;
        end
      else
        cond=0;
      end
    else
      cond=0;
    end
        
  end % MAIN 'while' cycle 
                   
    
  % OUTPUT VARIABLES FOR DATA(JJ,:)
    
  % Sobs1 (pDNA observed permeation)
  if condout == 1
    Sobs1(jj) = NaN;
    fprintf(['\nWarning:\nFunction did not find a\n'...
             'solution for row %s of input cell.\n'...
             'Please recheck the input values.\n'...
             'Also, Sobs may be too close to 1.\n']...
             ,strcounter)
  else
    Sobs1(jj) = fi1 * Cm1 / Cb1;
  end
    
  % Sobs2 (RNA observed permeation)
  if nargout >= 2
    if condout == 1
      Sobs2(jj) = NaN;
    else
      Sobs2(jj) = fi2*Cm2/Cb2;
    end
  end
    
  % Structure pro
  if nargout >= 3
    if condout == 1
      pro(len).dist      = 'No solution found';
      pro(len).pdna      = 'No solution found';
      pro(len).rna       = 'No solution found';
      pro(len).anion     = 'No solution found';
      pro(len).cation    = 'No solution found';
      pro(len).potential = 'No solution found';
    else
      pro(jj).dist      = y;
      pro(jj).pdna      = f(:, 1);
      pro(jj).rna       = f(:, 2);
      pro(jj).anion     = f(:, 3);
      pro(jj).cation    = f(:, 4);
      pro(jj).potential = f(:, 5);
    end
  end
    
  % Cell logF
  if nargout >= 4
    if condout == 1
      logF{jj+1, 1} = 'No solution found';
      logF{jj+1, 2} = 'No solution found';
      logF{jj+1, 3} = 'No solution found';
      logF{jj+1, 4} = 'No solution found';
      logF{jj+1, 5} = 'No solution found';
      logF{jj+1, 6} = 'No solution found';
      logF{jj+1, 7} = 'No solution found';
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
      logF{jj+1, 4} = Cm2;
      logF{jj+1, 5} = f(end, 2);
      logF{jj+1, 6} = Cm3;
      logF{jj+1, 7} = f(end, 3);
    end
  end
    
  % Plot polarization profiles
  % New figure for each row in data
  if nargout >= 5
    if condout == 1
      figstr = ['No solution found for values in row ',...
                num2str(counter)];
      figure('name',figstr);
      graph = true;
    else
      figstr = ['Polarization for values in row '...
                ,num2str(counter),' of input data'];
      figure('name',figstr)
      subplot(3, 2, 1)
      plot(y,f(:, 1), 'k-')
      title('pDNA')
      xlabel('y [m]')
      ylabel('C_1 [mol/m^3]')
      subplot(3, 2, 2)
      plot(y,f(:, 2),'k-')
      title(strRNA)
      xlabel('y [m]')
      ylabel('C_2 [mol/m^3]')
      subplot(3, 2, 3)
      plot(y,f(:, 3), 'k-')
      title(anion)
      xlabel('y [m]')
      ylabel('C_3 [mol/m^3]')
      subplot(3, 2, 4)
      plot(y,f(:, 4), 'k-')
      title(cation)
      xlabel('y [m]')
      ylabel('C_4 [mol/m^3]')
      subplot(3, 2, 5)
      plot(y,f(:, 5), 'k-')
      title('potential')
      xlabel('y [m]')
      ylabel('\Psi [V]')
      graph = true;
    end
  end
    
end % End of main 'for' cycle

strtoc = '\nElapsed time: %.6f seconds\n';
fprintf(strtoc, toc)