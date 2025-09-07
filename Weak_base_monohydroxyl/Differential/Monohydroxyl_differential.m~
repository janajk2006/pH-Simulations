% Obtained form the webpage: http://www.chemicalforums.com/index.php?topic=80995.0
tic
clear
clc

% Equilibrium constants
Kw = power(10,-14);
K1 = power(10,-4.6); % pka=14-pKb, 14-9.245=4.6


% Initial concentration of acid
C_acid = 1; % mol/dm^3
C_base = 1; % mol/dm^3

% Initial volume of acid
Acid_vol = 50; % ml
Base_vol = 10; % ml



% Graph settings
Titration_resolution=50; % number of titration points
dx = Acid_vol/(Titration_resolution-1);
x = 0.0:dx:Acid_vol;


% Array for all the species
BOH=size(Titration_resolution,1);
B=size(Titration_resolution,1);
H=size(Titration_resolution,1);
OH=size(Titration_resolution,1);
PH=size(Titration_resolution,1);
FBOH=size(Titration_resolution,1);
FHCl=size(Titration_resolution,1);



K=size(5);
Y_init=size(5);
Y_tol=size(5);

for i=1:5
Y_tol(i)=power(10,-16);
end

for i = 1:Titration_resolution

	% The total volume
	vol_tot = Base_vol + x(i);

	% Concentrations of acid and base
	C_a0 = (C_acid * x(i)) / vol_tot;
	C_b0 = (C_base * Base_vol) / vol_tot ;

	% Compute the concentration of hydrogen ions

		Y_init(1)=C_b0;
		Y_init(2)=0;
		Y_init(3)=C_a0;
		Y_init(4)=power(10,-7);
		Y_init(5)=power(10,-7);

		options = odeset('RelTol',1e-6,'AbsTol',Y_tol);



		K(1)=K1;
		K(2)=power(10,0);
		K(3)=power(10,1);			%6.3*power(10,-1);
		K(4)=power(10,-14);
		K(5)=power(10,0);


		[T,Y] = ode15s(@Monohydroxyl_ODE,[0,10000000],Y_init,options,K);


	S=size(Y);




	% Calculate all other species
	BOH(i,1) = Y(S(1),1);
	B(i,1) = Y(S(1),2);
	HCl(i,1)=Y(S(1),3);
	H(i,1) = Y(S(1),4);
	OH(i,1) = Y(S(1),5);
	FHCl(i,1)=C_a0;
	FBOH(i,1)=C_b0;

	% Calculate the pH value
	PH(i,1) = -log10(H(i,1));	

end

% Plot it
plot(x,PH,'red')
xlabel('Volume Acid (ml)')
ylabel('pH')
title('Titration (Volume of acid vs pH)')
toc
