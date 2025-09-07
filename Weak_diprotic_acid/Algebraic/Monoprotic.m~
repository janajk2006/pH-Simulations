% Obtained form the webpage: http://www.chemicalforums.com/index.php?topic=80995.0
clear;
clc;

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
for i = 1:Titration_resolution

	% The total volume
	vol_tot = Acid_vol + x(i);

	% Concentrations of acid and base
	C_a0 = (C_acid * Acid_vol) / vol_tot;
	C_b0 = (C_base * x(i)) / vol_tot ;

	% Compute the concentration of hydrogen ions
	% Solve the polynomial equation
	R = roots([1, (K1+C_b0), (C_b0*K1 - C_a0*K1 - Kw), -K1*Kw]);
	% Remove all chemically irrelevant solutions
	R = R(R==real(R) & R>0 & R<C_a0);
	% If there are more than one chemically relevant solution, display the number
	if not(length(R) == 1)
	disp(length(R))
	end
	H(i,1) = min(R(R==real(R) & R>0 & R<C_a0));
	% Calculate the pH value
	PH(i,1) = -log10(H(i,1));

	% Calculate all other species
	HA(i,1) = H(i,1)*C_a0/(K1+H(i,1));
	A(i,1) = K1*C_a0/(K1+H(i,1));
	OH(i,1) = Kw/H(i);
	FHA(i,1)=C_a0;
	FNa(i,1)=C_b0;
		

end

% Plot it
plot(x,PH,'red')
xlabel('Volume Alkali (ml)')
ylabel('pH')
title('Titration (Volume of alkali vs pH')
