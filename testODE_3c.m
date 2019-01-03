function [Cm1eOUT, Cm2eOUT, condout] = testODE_3c(Cm1, Cm2, fi, ...
    delta, x0, z1, z2, z3, F, R, T, Jv, D1, D2, D3)

Cm1IN=Cm1;
Cm2IN=Cm2;

condout=0;

Cp1 = fi * Cm1; Cp2 = Cm2; Cp3 = - (z1 * Cp1 + z2 * Cp2) / z3;
[~, f] = rk4arg(@nernst3c, [0 delta], x0, 100, Cp1, Cp2, Cp3, ...
  z1, z2, z3, F, R, T, Jv, D1, D2, D3);

res = iscres(f(:, 1));
if res == 1
  Cm1eOUT = Cm1IN; Cm2eOUT = Cm2IN;
  return
end

while res==0
  Cp1 = fi * Cm1; Cp2 = Cm2; Cp3 = - (z1 * Cp1 + z2 * Cp2) / z3;
  [~, f] = rk4arg(@nernst3c, [0 delta], x0, 100, Cp1, Cp2, ...
    Cp3, z1, z2, z3, F, R, T, Jv, D1, D2, D3);
  res = iscres(f(:, 1));
  dif2 = Cm1 - f(end, 1);
  if res == 1, break, end
  [mnt, expnt] = mantexpnt(Cm1);
  Cm1 = (mnt - 0.1) * 10 ^ (expnt);
end

if sign(dif2) == 1
  Cm1eOUT = Cm1; Cm2eOUT = Cm2;
  return
end

while sign(dif2) == -1
  [mnt, expnt] = mantexpnt(Cm1);
  Cm1 = (mnt + 0.01) * 10 ^ (expnt);
  Cp1 = fi * Cm1; Cp2 = Cm2; 
  Cp3 = - (z1 * Cp1 + z2 * Cp2) / z3;
  [~, f] = rk4arg(@nernst3c, [0 delta], x0, 100, Cp1, Cp2, ...
    Cp3, z1, z2, z3, F, R, T, Jv, D1, D2, D3);
  res = iscres(f(:, 1));
  if res == 0
    condout = 1;
    Cm1eOUT = Cm1; Cm2eOUT = Cm2;
    return
  end
  dif2 = Cm1 - f(end, 1);
end

Cm1eOUT = Cm1; Cm2eOUT = Cm2;