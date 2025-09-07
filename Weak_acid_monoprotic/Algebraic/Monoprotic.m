% Obtained form the webpage: http://www.chemicalforums.com/index.php?topic=80995.0
clear;
clc;


Max_species=1;
Max_roots=Max_species+2;
% Equilibrium constants
Kw = power(10,-14);
K1 = power(10,-4); % power(10,-4.76);


% Initial concentration of acid
C_acid = 1; % mol/dm^3
C_base = 1; % mol/dm^3

% Initial volume of acid
Acid_vol = 10; % ml
Base_vol = 25; % ml



% Graph settings
Titration_resolution=500; % number of titration points
dx = Base_vol/(Titration_resolution-1);
x = 0.0:dx:Base_vol;


% Array for all the species
HA=size(Titration_resolution,1);
A=size(Titration_resolution,1);
H=size(Titration_resolution,1);
PH=size(Titration_resolution,1);
OH=size(Titration_resolution,1);
FNa=size(Titration_resolution,1);
FHA=size(Titration_resolution,1);
Beta=size(Titration_resolution-1,1);

outfile=fopen('Results','w');

for j = 1:Max_roots
fprintf(outfile,'Roots(%d) ',j);
end
fprintf(outfile,'Acid_vol Base_vol Vol_tot HA(i,1) A(i,1) H(i,1) OH(i,1) PH(i,1) FHA(i,1) FNa(i,1) Beta(i,1)\n');


for i = 1:Titration_resolution

	% The total volume
	vol_tot = Acid_vol + x(i);

	% Concentrations of acid and base
	C_a0 = (C_acid * Acid_vol) / vol_tot;
	C_b0 = (C_base * x(i)) / vol_tot ;


	% Compute the concentration of hydrogen ions
	% Solve the polynomial equation

	% H^3 + (C_b0 + K1)*H^2 + (C_b0*K1 - C_a0*K1 - Kw)*H - K1*Kw

	R = roots([1, (C_b0 + K1), (C_b0*K1 - C_a0*K1 - Kw), - K1*Kw]);

	% Remove all chemically irrelevant solutions	
	check_index=1;
	Temp_H=0;
	for j = 1:Max_roots

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


	% Calculate the pH value
	PH(i,1) = -log10(H(i,1));

	% Calculate all other species
	%HA=(C_a0*H)/(H + K1)
	%A=(C_a0*K1)/(H + K1)

	HA(i,1) = (C_a0*H(i,1))/(H(i,1) + K1);
	A(i,1) = (C_a0*K1)/(H(i,1) + K1);
	OH(i,1) = Kw/H(i);
	FHA(i,1)=C_a0;
	FNa(i,1)=C_b0;
	if(i>1)
	Beta(i,1)=(FNa(i,1)-FNa(i-1,1))/(PH(i,1)-PH(i-1,1));
	else 
	Beta(i,1)=0;
	end
	
	fprintf(outfile,'%.12f %.12f %.12f %.12f %.12f %.20f %.20f %.12f %.12f %.12f %.12f\n',Acid_vol, x(i), vol_tot, HA(i,1),A(i,1),H(i,1), OH(i,1), PH(i,1),FHA(i,1),FNa(i,1),Beta(i,1));

end
fclose(outfile);

% Plot it
figure();
plot(x,PH,'red')
xlabel('Volume Alkali (ml)')
ylabel('pH')
title('Titration (Volume of alkali vs pH)')
