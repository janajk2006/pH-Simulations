% Obtained form the webpage: http://www.chemicalforums.com/index.php?topic=80995.0
tic
clear
clc

% Equilibrium constants
Kw = power(10,-14);
K1 = power(10,-4.76);


% Initial concentration of acid
C_acid = 1; % mol/dm^3
C_base = 1; % mol/dm^3

% Initial volume of acid
Acid_vol = 10; % ml
Base_vol = 50; % ml



% Graph settings
Titration_resolution=50; % number of titration points
dx = 50/(Titration_resolution-1);
x = 0.0:dx:Base_vol;


% Array for all the species
HA=size(Titration_resolution,1);
A=size(Titration_resolution,1);
H=size(Titration_resolution,1);
PH=size(Titration_resolution,1);
OH=size(Titration_resolution,1);
FNa=size(Titration_resolution,1);
FHA=size(Titration_resolution,1);


K=size(6);
Y_init=size(4);
Y_tol=power(10,-4);
for i = 1:Titration_resolution

	% The total volume
	vol_tot = Acid_vol + x(i);

	% Concentrations of acid and base
	C_a0 = (C_acid * Acid_vol) / vol_tot;
	C_b0 = (C_base * x(i)) / vol_tot ;

		Y_init(1)=C_a0;
		Y_init(2)=0;
		Y_init(3)=1e-7;
		Y_init(4)=1e-3;



		K(1)=power(10,-5);
		K(2)=power(10,0);
		K(3)=power(10,-14);			%6.3*power(10,-1);
		K(4)=power(10,0);
		K(5)=C_a0;
		K(6)=C_b0;

	% Compute the concentration of hydrogen ions
M = [1 0 0 0
     0 1 0 0
     0 0 1 0
     0 0 0 0];


		options = odeset('Mass',M,'RelTol',1e-4, 'AbsTol',[1e-16 1e-16  1e-16  1e-16]);

		[T,Y] = ode15s(@Monoprotic_DAE_func,[0 10000000],Y_init,options,K);


	S=size(Y);




	% Calculate all other species
	HA(i,1) = Y(S(1),1);
	A(i,1) = Y(S(1),2);
	H(i,1) = Y(S(1),3);
	%OH(i,1) = Y(S(1),5);
	FHA(i,1)=C_a0;
	FNa(i,1)=C_b0;

	% Calculate the pH value
	PH(i,1) = -log10(H(i,1));	

end

% Plot it
plot(x,PH,'red')
xlabel('Volume Alkali (ml)')
ylabel('pH')
title('Titration (Volume of alkali vs pH')
toc
