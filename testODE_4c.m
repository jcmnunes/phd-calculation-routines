function [Cm1eOUT, Cm2eOUT, Cm3eOUT, condout] = testODE_4c(Cm1, Cm2, Cm3,...
    fi1, fi2, delta, x0, z1, z2, z3, z4, F, R, T, Jv, D1, D2, ...
    D3, D4, strcounter)

Cm1IN = Cm1;
Cm2IN = Cm2;
Cm3IN = Cm3;

condout = 0;

Cp1 = fi1 * Cm1; Cp2 = fi2 * Cm2; Cp3 = Cm3; 
Cp4 = - (z1 * Cp1 + z2 * Cp2 + z3 * Cm3) / z4;
[~, f] = rk4arg(@nernst4c, [0 delta], x0, 100, Cp1, Cp2, ...
  Cp3, Cp4, z1, z2, z3, z4, F, R, T, Jv, D1, D2, D3, D4);

res = iscres(f(:, 1));
if res == 1
  Cm1eOUT = Cm1IN; Cm2eOUT = Cm2IN; Cm3eOUT = Cm3IN;
  return
end

ct = 1;
while res == 0
  Cp1 = fi1 * Cm1; Cp2 = fi2 * Cm2; Cp3 = Cm3;
  Cp4 = -(z1 * Cp1 + z2 * Cp2 + z3 * Cp3) / z4;
  [~, f] = rk4arg(@nernst4c, [0 delta], x0, 100, Cp1, Cp2, ...
    Cp3, Cp4, z1, z2, z3, z4, F, R, T, Jv, D1, D2, D3, D4);
  res = iscres(f(:, 1));
  dif2 = Cm1 - f(end, 1);
  if res == 1, break, end
  [mnt, expnt] = mantexpnt(Cm1);
  Cm1 = (mnt - 0.1) * 10 ^ (expnt);
  % Prevent infinite loop
  if ct == 500
    fprintf(['\nWarning:\nFunction is taking too long\n'...
             'to find a solution for row %s of input cell.\n'...
             'Use ctrl + c to abort and try a lower flux value.\n']...
             ,strcounter)
  end
  if ct == 501
    fprintf('\n*** The function may abort soon ***\n');
  end
  if ct == 1500 
    condout = 1;
    fprintf('\n*** Aborting ... ***\n')   
  end
  if ct == 2000 
    condout = 1;
    Cm1eOUT = Cm1; Cm2eOUT = Cm2; Cm3eOUT = Cm3;
    return 
  end 
  ct = ct + 1;
end

if sign(dif2) == 1
  Cm1eOUT=Cm1;Cm2eOUT=Cm2;Cm3eOUT=Cm3;
  return
end

while sign(dif2) == -1
  [mnt, expnt] = mantexpnt(Cm1);
  Cm1 = (mnt + 0.01) * 10 ^ (expnt);
  Cp1 = fi1 * Cm1; Cp2 = fi2 * Cm2; Cp3 = Cm3;
  Cp4 = - (z1 * Cp1 + z2 * Cp2 + z3 * Cp3) / z4;
  [~, f] = rk4arg(@nernst4c, [0 delta], x0, 100, Cp1, Cp2, ...
    Cp3, Cp4, z1, z2, z3, z4, F, R, T, Jv, D1, D2, D3, D4);
  res = iscres(f(:, 1));
  if res == 0
    condout = 1;
    Cm1eOUT = Cm1; Cm2eOUT = Cm2; Cm3eOUT = Cm3;
    return
  end
  dif2 = Cm1 - f(end, 1);
end

Cm1eOUT = Cm1; Cm2eOUT = Cm2; Cm3eOUT = Cm3;