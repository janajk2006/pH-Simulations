% Obtained form the webpage: http://www.chemicalforums.com/index.php?topic=80995.0
tic
clear
clc

% Equilibrium constants
Kw = power(10,-14);
K1 = power(10,-8.1); % pka=14-pKb, 14-9.245=4.6


% Initial concentration of acid
C_acid = .5; % mol/dm^3
C_base = .49997; % mol/dm^3

% Initial volume of acid
Acid_vol = 1; % ml
Base_vol = 2; % ml



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
R =size(3,1);

outfile=fopen('Results','w');
fprintf(outfile,'Roots(3) BOH(i,1) B(i,1) H(i,1) OH(i,1) PH(i,1) FHCl(i,1) FBOH(i,1)\n');

for i = 1:Titration_resolution

	% The total volume
	vol_tot = Base_vol + x(i);

	% Concentrations of acid and base
	C_a0 = (C_acid * x(i)) / vol_tot;
	C_b0 = (C_base * Base_vol) / vol_tot ;

	% Compute the concentration of hydrogen ions
	% Solve the polynomial equation
	R = roots([K1, (K1*C_b0+Kw-C_a0*K1), -(K1*Kw + C_a0*Kw), -Kw*Kw]);

	% Remove all chemically irrelevant solutions	
	check_index=1;
	Temp_H=0;
	for j = 1:3
		fprintf(outfile,'%.20f ',R(j,1));
		if R(j,1)==real(R(j,1)) & R(j,1)>0 & R(j,1)<C_a0
			if check_index==1
			Temp_H=R(j,1);
			check_index=2;
			end

			if Temp_H >R(j,1) 
			Temp_H=R(j,1);	
			end
		end
	end

			H(i,1)=Temp_H;	




	% Calculate all other species
	BOH(i,1) = Kw*C_b0/(K1*H(i,1)+Kw);
	B(i,1) = K1*H(i,1)*C_b0/(K1*C_a0+Kw);
	HCl(i,1)=C_a0;
	OH(i,1) = Kw/H(i,1);
	FHCl(i,1)=C_a0;
	FBOH(i,1)=C_b0;

	% Calculate the pH value
	PH(i,1) = -log10(H(i,1));


	fprintf(outfile,'%.12f %.12f %.20f %.20f %.12f %.12f %.12f\n',BOH(i,1),B(i,1),H(i,1), OH(i,1), PH(i,1),FHCl(i,1),FBOH(i,1));
end
fclose(outfile);

% Plot it
plot(x,PH,'red')
xlabel('Volume Acid (ml)')
ylabel('pH')
title('Titration (Volume of acid vs pH)')
toc
