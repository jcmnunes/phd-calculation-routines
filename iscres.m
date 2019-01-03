function res = iscres(F)
% RES = ISCRES(F) returns 0 if test fails
% and 1 otherwise. Called in functions np_pdna_3c,
% np_rna_3c and np_4c.
%
% see also np_pdna_3c, np_rna_3c, np_4c

% Author: J. Nunes <josenunes34@gmail.com>
% Date: 2013

res = 1;
d = diff(F);
if any(d < 0)
  d(d > 0) = 1;
  d(d < 0) = 0;
  sign_changes = find(diff(d) ~= 0);
  len = length(sign_changes);
  smooth = find(abs(diff(F)) > 2 * mean(abs(diff(F))), 1);
  if len == 1
    if isempty(smooth) 
      res = 1;
      return
    else
      res = 0;
    end
  else
   res = 0;
  end
end




















