
% Decodify bitstrings
function x = decode8(v,M);
% x		-> real value (precision: 6)
% v		-> binary string (length: 22)
% M     -> Numero maximo que deseo obtener,en rango creciente de 1 hasta M.
v = fliplr(v); %Invierte la matriz que recibe. 
s = size(v);%Calcula el tamaño de la matriz que recibe.
aux = 0:1:7; %Crea un auxiliar de 22 bits.
aux = ones(s(1),1)*aux; %Multiplica un vector columna de 100x1, por el vector aux de 1x22.
x1 = sum((v.*2.^aux)');%Multiplica cada uno de los elementos de la matriz que recibe por 2, 
%ese resultado se eleva al valor correspondiente que haya en aux;al final,
%transpone todo el resultado,suma cada renglon para dejar un vector 1x100.
x =  round(x1 .* (M / 255));%Hago la proporcion del numero maximo real que deseo obtener 
%entre el numero maximo que puedo representar con 22 bits(4194303),y la
%multiplico elemento a elemento por el vector anterior 1x100,y redondeo.