%Declara��o de vari�veis
% m %Numero de esta��es na envergadura nas quais se encontram a circula��o e downwash
% n %Integrador que define a esta��o na linha de quarto de corda da asa para o qual � calculado o valor de circula��o
% v %Integrador que define um ponto especifico dentro da forma em planta da asa para o qual a condi��o de fronteira que especifica o n�o atravessamento de escoamento na asa � aplicada
% eta %Coordenada lateral adimensional medida perpendicular ao plano de simetria, frac��o da semienvergadura.
% eta_cp %Localiza��o do centro de press�o na envergadura
% eta_cpa %Localiza��o do centro de press�o na envergadura devido � sustenta��o adicional na asa
% fi=acos(eta); %Coordenada trigonom�trica na envergadura (radianos)
% fi_n=acos(n*pi)/8; %Valor de fi na esta��o n
% alfa %�ngulo de ataque da asa medido no plano paralelo ao plano de simetria (graus ou radianos)
% alfa_v %�ngulo de ataque na esta��o v (graus ou radianos)
% alfa_v0 %Angulo de ataque na esta��o v para sustenta��o neta nula na asa (graus ou radianos)
% alfa_r0 %Angulo de ataque da raiz da asa para sustenta��o neta nula na asa (graus ou radianos)
% epsilon_v=alfa_v0-alfa_r0; %�ngulo de torcimento da asa na esta��o v medida relativa � se��o da raiz da asa
% GAMA %Circula��o (p�s quadrados por segundo)
% b %Envergadura da asa medida perpendicular do plano de simetria (p�s)
% V %Velocidade do escoamento livre (p�s por segundo)
% Gn=(GAMA/b*V); %Circula��o adimensional (GAMA/b*V), id�ntica ao coeficiente de carga cl*c/2*b, na esta��o da envergadura n
% Gn0 %Valor de Gn para sustenta��o neta nula na asa
% Gna %Valor de Gn devido � sustenta��o adicional na asa
% A_v_n %Coeficiente que depende da geometria da asa e que indica a influencia do carregamento arbitr�rio na esta��o n da envergadura no angulo de downwash na esta��o da envergadura v
% a_v_n %Coeficiente que depende da geometria da asa e que indica a influencia da carga sim�trica na esta��o n no angulo de downwash na esta��o v
% c %Corda local da asa medida paralela ao plano de simetria (p�s)
% c_v %Corda local no ponto v (p�s)
% c_ag=S/b; %Corda prom�dio
% MAC %Corda media aerodin�mica (p�s)
% lambda %Conicidade
% A=(b^2)/S; %Afilamento
% S %Superf�cie da asa (p�s quadrados)
% LAMBDA %�ngulo de enflechamento geom�trico a um quarto de corda, positivo quando o enflechamento � para traz da linha normal do plano de simetria (graus)
% M %N�mero de Mach
% BETA=sqrt((1-M^2)) %Par�metro de compressibilidade 
% delta_beta=atan((tan(LAMBDA))/BETA); %par�metro de enflechamento compress�vel 
% cl %Coeficiente de sustenta��o do perfil
% cl_alfa %Inclina��o da curva de sustenta��o do perfil (radianos ou graus)
% cl_0 %Coeficiente de sustenta��o do perfil para sustenta��o neta nula na asa
% cl_alfa %Coeficiente de sustenta��o do perfil devido � sustenta��o adicional na asa
% CL %coeficiente de sustenta��o da asa
% CL_alfa %inclina��o da curva de sustenta��o (radianos ou graus)
% Cm %Coeficiente de momento 
% Cm_0 %Coeficiente de momento para sustenta��o neta nula da asa
% Cm %Coeficiente de momento devido � sustenta��o adicional na asa
% CD_i %Coeficiente de arrastro induzido
% CD_i0 %Coeficiente de arrasto induzido devido � carga b�sica (sustenta��o neta nula na asa)
% CD_i_alfa %Coeficiente de arrastro induzido devido � carga adicional (sustenta��o neta na asa)
% a_c %Posi��o longitudinal do centro aerodin�mico, medido do borde de ataque da MAC em porcentagem da MAC
% GAMA %Circula��o (p�s quadrados por segundo)
% h %Distancia absoluta de um v�rtice at� um ponto de downwash, medido perpendicular ao v�rtice (p�s)
% rho %densidade do ar (slugs por p� cubico)
% w %velocidade induzida, perpendicular � linha da corda media da asa, positiva para downwash (p� por segundo)
% q %Press�o din�mica (libra por p� quadrado)
% k %rela��o da inclina��o da curva de sustenta��o experimental do perfil e a valor te�rico de 2*pi/BETA, ambos tomados ao mesmo n�mero de Mach
% k_v %Valor de k na esta��o v
% d_v %fator de escala
% H_v=d_v*(1/k_v)*(b/(c_v/BETA)) %Par�metro geom�trico da asa d_v*(1/k_v)*(b/(c_v/BETA))

n=[1 2 3 4 5 6 7 8];
eta=cos((n*pi)/8)

