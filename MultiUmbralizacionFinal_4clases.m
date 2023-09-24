%%%%Artificial Inmune Systems-Aplicacion multiumbralizacion imagen a
%%%%escala de grises. Version con tres gaussianas.
%%%%Ingeniero Jose Valentin Osuna Enciso, Marzo,Abril,2009
%%%%%Maestria en Ciencias en Ingenieria, CUCEI.
%%%%%%Version final con 3 gaussianas.
clear all
format long
%Carga imagen,convierte a escala de grises,saca su histograma,calcula maximo.
DB=imread('243095.jpg');
% I2.jpg ->3 clases, quizas 4.
% I3.jpg ->2 clases.
% I4.jpg ->2 clases.
% I5.jpg ->3 clases, quizas 4.
% I6.jpg ->3 clases, quizas 4 o 5.
% I7.jpg ->3 clases.
% I8.jpg ->2 clases, quizas 3
% I9.jpg ->2 clases, quizas 3 o 4.
% I10.jpg ->2 clases, quizas 3 o 4.
% I11.tif ->3 clases
numClases=4;DB=DB(:,:,1);%DB=rgb2gray(DB);%Este rgb2gray se elimina si la imagen es tif.
DB4=DB;
h=imhist(DB);
[m3 a3]=max(h);

%%%%%%%%%%Aqui empieza el primer CLONALG que busca las medias%%%%%%%%%%%%
V=cadeia(200,16,0,0,0);
x=0:255; n = size(V,1); it = 0;
x1=decode8(V(:,1:8),255);x2=decode8(V(:,9:16),255);
for ind1=1:n
    if(x2(ind1)>=x1(ind1))
        [fit1(ind1) ubica(ind1)]=max(h(x1(ind1):x2(ind1),1));
        ubica(ind1)=ubica(ind1)+x1(ind1)-1;
    else
        [fit1(ind1) ubica(ind1)]=max(h(x2(ind1):x1(ind1),1));
        ubica(ind1)=ubica(ind1)+x2(ind1)-1;
    end
end
[a ind]=sort(fit1);
ubicaOrden=ubica(ind(end-n+1:end));
for ind1=1:n
    temp=ubicaOrden(1,ind1);
    ubicaOrden(2,ind1)=1;
    for ind2=ind1+1:n
        if(temp==ubicaOrden(1,ind2)&& ubicaOrden(1,ind2)>0)
            ubicaOrden(1,ind2)=0;
            ubicaOrden(2,ind1)=ubicaOrden(2,ind1)+1;
        end
    end
end 
[a ind]=sort(ubicaOrden(2,:));         %Nota:uso las mismas variables para no crear otras.
ubicaOrden=ubicaOrden(1,ind(end-n+1:end));%igual aqui.Los mejores quedan al final de ubicaOrden.
Med=ubicaOrden(1,end-numClases+1:end);Med=sort(Med);A=h(Med(1,1:numClases));
valle=min(h((Med(1,1)):(Med(1,2))));
valle2=find(h(:,1)==valle); %En 'M' se guardan las medias, en 'A' las amplitudes.
%Aqui termina el primer CLONALG que busca las medias y amplitudes maximas.%%%%%%%%%%%%%%%

%%%%%%%%%%%%%Aqui empieza el segundo CLONALG que busca las desviaciones estandar:%%%%%%
V=cadeia(60,8*8,0,0,0); caso=1;salida=1;
x=0:255; gen = 200; n = size(V,1); fat = .2; [N,L] = size(V); it = 0; 
while it<gen && salida<10
disp(sprintf('Numero de generacion = %d, salida= %d, Fitness= %d',it, salida,caso));  
desvest1 = decode8(V(:,1:8),30); desvest2 = decode8(V(:,9:16),30); desvest3 = decode8(V(:,17:24),30);desvest4 = decode8(V(:,25:32),30);
Ampli1=decode8(V(:,33:40),A(1,1)); Ampli2 = decode8(V(:,41:48),A(2,1)); Ampli3 = decode8(V(:,49:56),A(3,1));Ampli4 = decode8(V(:,57:64),A(4,1));
T = []; cs = [];
%
for ind1=1:n
   gaussianas(:,ind1)=(Ampli1(1,ind1)*exp(-((x-Med(1,1)).^2)/(2*(desvest1(ind1)^2))))+(Ampli2(1,ind1)*exp(-((x-Med(1,2)).^2)/(2*(desvest2(ind1)^2))))+(Ampli3(1,ind1))*exp(-((x-Med(1,3)).^2)/(2*(desvest3(ind1)^2)))+(Ampli4(1,ind1))*exp(-((x-Med(1,4)).^2)/(2*(desvest4(ind1)^2)));
   errortotal(:,ind1)=((h-gaussianas(:,ind1))).^2;
   error(ind1)=sqrt(sum(errortotal(:,ind1))/256);
end

[a ind]=sort(error);
if(caso==(a(1)))
    salida=salida+1;
else
    salida=0;
end
%%%%%%%%%%%%Termino el calculo de los errores:%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
caso=a(1);
valDE1 = desvest1(ind(end-n+1:end));valDE2 = desvest2(ind(end-n+1:end)); valDE3 = desvest3(ind(end-n+1:end));valDE4 = desvest4(ind(end-n+1:end));
valA1=Ampli1(ind(end-n+1:end));valA2 = Ampli2(ind(end-n+1:end)); valA3 = Ampli3(ind(end-n+1:end));valA4 = Ampli4(ind(end-n+1:end));
% Reproduction
[T,pcs] = reprod(n,fat,N,ind,V,T);
% Hypermutation
M = rand(size(T,1),L) <= 0.35;
T = T - 2 .* (T.*M) + M;
T(pcs,:) = V(fliplr(ind(end-n+1:end)),:);%En las posiciones 'pcs' de T, coloca los individuos de C basados en afinidad.
%%%%%%%Busco los valores minimos:%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
desvest1 = decode8(T(:,1:8),30); desvest2 = decode8(T(:,9:16),30); desvest3 = decode8(T(:,17:24),30);desvest4 = decode8(T(:,25:32),30);pcs = [0 pcs];
Ampli1=decode8(T(:,33:40),A(1,1)); Ampli2 = decode8(T(:,41:48),A(2,1)); Ampli3 = decode8(T(:,49:56),A(3,1));Ampli4 = decode8(T(:,57:64),A(4,1));
for ind1=1:size(T,1)
   gaussianas2(:,ind1)=(Ampli1(1,ind1)*exp(-((x-Med(1,1)).^2)/(2*(desvest1(ind1)^2))))+(Ampli2(1,ind1)*exp(-((x-Med(1,2)).^2)/(2*(desvest2(ind1)^2))))+(Ampli3(1,ind1))*exp(-((x-Med(1,3)).^2)/(2*(desvest3(ind1)^2)))+(Ampli4(1,ind1))*exp(-((x-Med(1,4)).^2)/(2*(desvest4(ind1)^2)));
   errortotal2(:,ind1)=((h-gaussianas2(:,ind1))).^2;
   fit2(ind1)=sqrt(sum(errortotal2(:,ind1))/256);
end
for ind1=1:n,
   [out(ind1),bcs(ind1)] = min(fit2(pcs(ind1)+1:pcs(ind1+1)));		% Problema minimizacion del error
   bcs(ind1) = bcs(ind1) + pcs(ind1);
   if ind1==(n-1)
       ind1=ind1;
   end
end;
 %%%%%% % %  %%%Mezcla de individuos.%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   V(fliplr(ind(end-n+1:end)),:) = T(bcs,:);
%%%%%%%%%%%Termina mezcla individuos.%%%%%%%%%%%%%%%%%%%%%%%%%%%    
it=it+1;
end;
%%%%%%Termina segundo CLONALG.%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
DE1=valDE1(1);DE2=valDE2(1);DE3=valDE3(1);DE4=valDE4(1);A1=valA1(1);A2=valA2(1);A3=valA3(1);A4=valA4(1);
Resultado=(A1*exp(-((x-Med(1,1)).^2)/(2*(DE1^2))))+(A2*exp(-((x-Med(1,2)).^2)/(2*(DE2^2))))+(A3*exp(-((x-Med(1,3)).^2)/(2*(DE3^2))))+(A4*exp(-((x-Med(1,4)).^2)/(2*(DE4^2))));
plot(h,'r')
hold on
plot(Resultado,'k'),figure,
plot((A4*exp(-((x-Med(1,4)).^2)/(2*(DE4^2)))),'g'),hold on
plot((A3*exp(-((x-Med(1,3)).^2)/(2*(DE3^2)))),'g')
plot((A2*exp(-((x-Med(1,2)).^2)/(2*(DE2^2)))),'g')
plot((A1*exp(-((x-Med(1,1)).^2)/(2*(DE1^2)))),'g')
plot(Resultado,'k')
%%%%%Realizo umbralizacion imagen escala de grises:%%%%%%%%%%%%%%%%%%%%%%%
a1=(DE1^2)-(DE2^2);
a2=(DE2^2)-(DE3^2);
a3=(DE3^2)-(DE4^2);
b1=2*((Med(1,1)*(DE2^2))-(Med(1,2)*(DE1^2)));
b2=2*((Med(1,2)*(DE3^2))-(Med(1,3)*(DE2^2)));
b3=2*((Med(1,3)*(DE4^2))-(Med(1,4)*(DE3^2)));
c1=((DE1*Med(1,2))^2)-((DE2*Med(1,1))^2)+(2*((DE1*DE2)^2)*log((DE2*A(1,1))/(DE1*A(2,1))));
c2=((DE2*Med(1,3))^2)-((DE3*Med(1,2))^2)+(2*((DE3*DE2)^2)*log((DE3*A(2,1))/(DE2*A(3,1))));
c3=((DE3*Med(1,4))^2)-((DE4*Med(1,3))^2)+(2*((DE4*DE3)^2)*log((DE4*A(3,1))/(DE3*A(4,1))));
T1a=(-b1+sqrt((b1^2)-(4*a1*c1)))/(2*a1);
T1b=(-b1-sqrt((b1^2)-(4*a1*c1)))/(2*a1);
T2a=(-b2+sqrt((b2^2)-(4*a2*c2)))/(2*a2);
T2b=(-b2-sqrt((b2^2)-(4*a2*c2)))/(2*a2);
T3a=(-b3+sqrt((b3^2)-(4*a3*c3)))/(2*a3);
T3b=(-b3-sqrt((b3^2)-(4*a3*c3)))/(2*a3);
[fila columna]=size(DB);
for ind1=1:fila
    for ind2=1:columna
        if (DB(ind1,ind2)<=T1b)&&(DB(ind1,ind2)>=0)
            DBsegmented(ind1,ind2)=0;
        elseif (DB(ind1,ind2)<=T2b)&&(DB(ind1,ind2)>T1b)
            DBsegmented(ind1,ind2)=0.25;
        elseif(DB(ind1,ind2)<=T3b)&&(DB(ind1,ind2)>T2b)
            DBsegmented(ind1,ind2)=0.5;
        elseif(DB(ind1,ind2)>T3b)
            DBsegmented(ind1,ind2)=1;
        end
    end
end
figure,DBsegmented=mat2gray(DBsegmented);
imshow(DBsegmented)
