

% Control of pm
function [pm] = pmcont(pm,pma,pmr,it,itpm);
% pma	-> initial value
% pmr	-> control rate
% itpm	-> iterations for restoring
if rem(it,itpm) == 0,
   pm = pm * pmr;
   if rem(it,10*itpm) == 0,
      pm = pma;
   end;
end;