% Function CADEIA
function [ab,ag] = cadeia(n1,s1,n2,s2,bip)
if nargin == 2,
   n2 = n1; s2 = s1; bip = 1;
elseif nargin == 4,
   bip = 1;
end;
% Antibody (Ab) chains.
ab = 2 .* rand(n1,s1) - 1;%Crea cadenas de anticuerpos,aleatorias,de n1xs1,las multiplica elemento a elemento por 2 y al final les resta 1 a cada elemento.
if bip == 1,
   ab = hardlims(ab);
else,
   ab = hardlim(ab);%Convierte la matriz n1xs1 a binaria.
end;
% Antigen (Ag) chains,antigenos generan la formacion de anticuerpos.
ag = 2 .* rand(n2,s2) - 1;
if bip == 1,
   ag = hardlims(ag);
else,
   ag = hardlim(ag);
end;
% End Function CADEIA