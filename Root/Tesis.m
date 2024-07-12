%Declaração de variáveis
% m %Numero de estações na envergadura nas quais se encontram a circulação e downwash
% n %Integrador que define a estação na linha de quarto de corda da asa para o qual é calculado o valor de circulação
% v %Integrador que define um ponto especifico dentro da forma em planta da asa para o qual a condição de fronteira que especifica o não atravessamento de escoamento na asa é aplicada
% eta %Coordenada lateral adimensional medida perpendicular ao plano de simetria, fracção da semienvergadura.
% eta_cp %Localização do centro de pressão na envergadura
% eta_cpa %Localização do centro de pressão na envergadura devido à sustentação adicional na asa
% fi=acos(eta); %Coordenada trigonométrica na envergadura (radianos)
% fi_n=acos(n*pi)/8; %Valor de fi na estação n
% alfa %Ângulo de ataque da asa medido no plano paralelo ao plano de simetria (graus ou radianos)
% alfa_v %Ângulo de ataque na estação v (graus ou radianos)
% alfa_v0 %Angulo de ataque na estação v para sustentação neta nula na asa (graus ou radianos)
% alfa_r0 %Angulo de ataque da raiz da asa para sustentação neta nula na asa (graus ou radianos)
% epsilon_v=alfa_v0-alfa_r0; %Ângulo de torcimento da asa na estação v medida relativa à seção da raiz da asa
% GAMA %Circulação (pês quadrados por segundo)
% b %Envergadura da asa medida perpendicular do plano de simetria (pês)
% V %Velocidade do escoamento livre (pês por segundo)
% Gn=(GAMA/b*V); %Circulação adimensional (GAMA/b*V), idêntica ao coeficiente de carga cl*c/2*b, na estação da envergadura n
% Gn0 %Valor de Gn para sustentação neta nula na asa
% Gna %Valor de Gn devido à sustentação adicional na asa
% A_v_n %Coeficiente que depende da geometria da asa e que indica a influencia do carregamento arbitrário na estação n da envergadura no angulo de downwash na estação da envergadura v
% a_v_n %Coeficiente que depende da geometria da asa e que indica a influencia da carga simétrica na estação n no angulo de downwash na estação v
% c %Corda local da asa medida paralela ao plano de simetria (pês)
% c_v %Corda local no ponto v (pês)
% c_ag=S/b; %Corda promédio
% MAC %Corda media aerodinâmica (pês)
% lambda %Conicidade
% A=(b^2)/S; %Afilamento
% S %Superfície da asa (pês quadrados)
% LAMBDA %Ângulo de enflechamento geométrico a um quarto de corda, positivo quando o enflechamento é para traz da linha normal do plano de simetria (graus)
% M %Número de Mach
% BETA=sqrt((1-M^2)) %Parâmetro de compressibilidade 
% delta_beta=atan((tan(LAMBDA))/BETA); %parâmetro de enflechamento compressível 
% cl %Coeficiente de sustentação do perfil
% cl_alfa %Inclinação da curva de sustentação do perfil (radianos ou graus)
% cl_0 %Coeficiente de sustentação do perfil para sustentação neta nula na asa
% cl_alfa %Coeficiente de sustentação do perfil devido à sustentação adicional na asa
% CL %coeficiente de sustentação da asa
% CL_alfa %inclinação da curva de sustentação (radianos ou graus)
% Cm %Coeficiente de momento 
% Cm_0 %Coeficiente de momento para sustentação neta nula da asa
% Cm %Coeficiente de momento devido à sustentação adicional na asa
% CD_i %Coeficiente de arrastro induzido
% CD_i0 %Coeficiente de arrasto induzido devido à carga básica (sustentação neta nula na asa)
% CD_i_alfa %Coeficiente de arrastro induzido devido à carga adicional (sustentação neta na asa)
% a_c %Posição longitudinal do centro aerodinâmico, medido do borde de ataque da MAC em porcentagem da MAC
% GAMA %Circulação (pês quadrados por segundo)
% h %Distancia absoluta de um vórtice até um ponto de downwash, medido perpendicular ao vórtice (pês)
% rho %densidade do ar (slugs por pé cubico)
% w %velocidade induzida, perpendicular à linha da corda media da asa, positiva para downwash (pé por segundo)
% q %Pressão dinâmica (libra por pé quadrado)
% k %relação da inclinação da curva de sustentação experimental do perfil e a valor teórico de 2*pi/BETA, ambos tomados ao mesmo número de Mach
% k_v %Valor de k na estação v
% d_v %fator de escala
% H_v=d_v*(1/k_v)*(b/(c_v/BETA)) %Parâmetro geométrico da asa d_v*(1/k_v)*(b/(c_v/BETA))

n=[1 2 3 4 5 6 7 8];
eta=cos((n*pi)/8)

