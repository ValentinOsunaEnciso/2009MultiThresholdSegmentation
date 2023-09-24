
% Reproduction
function [T,pcs] = reprod(n,fat,N,ind,P,T);
% n		-> number of clones
% fat	-> multiplying factor
% ind	-> best individuals
% T		-> temporary population
% pcs	-> final position of each clone
if n == 1,
   cs = N;
   T = ones(N,1) * P(ind(1),:);
else,
   for i=1:n,
      % cs(i) = round(fat*N/i);
      cs(i) = round(fat*N);
      pcs(i) = sum(cs);
      T = [T; ones(cs(i),1) * P(ind(end-i+1),:)];
   end;
end;
