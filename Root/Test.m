
AR=8;
con=1;
% enfle=-45;
m=201;

Amin=-50;  %Angulo minimo (graus).
Amax=50;  %Angulo max (graus).
Adel=5;   %Incremento.
A = Amin:Adel:Amax;
for i=1:length(A);
     alfaA = A(i);
[Gma,fi_n]=slope(m,alfaA,AR,con);
incli(i)= ((pi*AR)/(m+1))*(Gma((m+1)/2,:)+sum(Gma(1:1:(m-1)/2).*(2*(sin(fi_n(1:1:(m-1)/2)))')));
end



plot(A,incli,'k')
axis([-50,50,0,7])
 hold on
grid  on
grid minor
 hold on
 xlabel('Enfle')
 ylabel('Slope')
