% Obtained form the webpage: http://www.chemicalforums.com/index.php?topic=80995.0
clear;
clc;

Max_species=3;
Max_roots=Max_species+2;
% Equilibrium constants
Kw = power(10,-14);
K1 = power(10,-3); % power(10,-4.76);
K2 = power(10,-6); % power(10,-7);
K3 = power(10,-9); % power(10,-7);

% Initial concentration of acid
C_acid = 1; % mol/dm^3
C_base = 1; % mol/dm^3

% Initial volume of acid
Acid_vol = 10; % ml
Base_vol = 45; % ml



% Graph settings
Titration_resolution=500; % number of titration points
dx = Base_vol/(Titration_resolution-1);
x = 0.0:dx:Base_vol;


% Array for all the species
H3A=size(Titration_resolution,1);
H2A=size(Titration_resolution,1);
HA=size(Titration_resolution,1);
A=size(Titration_resolution,1);
H=size(Titration_resolution,1);
PH=size(Titration_resolution,1);
OH=size(Titration_resolution,1);
FNa=size(Titration_resolution,1);
FH3A=size(Titration_resolution,1);

outfile=fopen('Results','w');

for j = 1:Max_roots
fprintf(outfile,'Roots(%d) ',j);
end
fprintf(outfile,'Acid_vol  x(i) vol_tot H3A(i,1) H2A(i,1) HA(i,1) A(i,1) H(i,1) OH(i,1) PH(i,1) FH3A(i,1) FNa(i,1) Beta(i,1)\n');


for i = 1:Titration_resolution

	% The total volume
	vol_tot = Acid_vol + x(i);

	% Concentrations of acid and base
	C_a0 = (C_acid * Acid_vol) / vol_tot;
	C_b0 = (C_base * x(i)) / vol_tot ;

	% Compute the concentration of hydrogen ions
	% Solve the polynomial equation
	%H^5 + (C_b0 + K1)*H^4 + (C_b0*K1 - C_a0*K1 - Kw + K1*K2)*H^3 + (C_b0*K1*K2 - 2*C_a0*K1*K2 - K1*Kw + K1*K2*K3)*H^2 + (C_b0*K1*K2*K3 - 3*C_a0*K1*K2*K3 - K1*K2*Kw)*H - K1*K2*K3*Kw =0;


	R = roots([1, (C_b0 + K1), (C_b0*K1 - C_a0*K1 - Kw + K1*K2), (C_b0*K1*K2 - 2*C_a0*K1*K2 - K1*Kw + K1*K2*K3), (C_b0*K1*K2*K3 - 3*C_a0*K1*K2*K3 - K1*K2*Kw), - K1*K2*K3*Kw]);
	% Remove all chemically irrelevant solutions	
	check_index=1;
	Temp_H=0;
	for j = 1:Max_roots
		fprintf(outfile,'%.20f ',R(j,1));
		if R(j,1)==real(R(j,1)) && R(j,1)>0 && R(j,1)<C_a0
			if check_index==1
			Temp_H=R(j,1);
			check_index=2;
			else 
				if Temp_H >= R(j,1) 
				Temp_H=R(j,1);	
				end
			end
		end
	end

			H(i,1)=Temp_H;	


	% Calculate the pH value
	PH(i,1) = -log10(H(i,1));

	% Calculate all other species
	%fA =(C_a0*K1*K2*K3)/(H^3 + K1*H^2 + K1*K2*H + K1*K2*K3)
	%fHA =(C_a0*H*K1*K2)/(H^3 + K1*H^2 + K1*K2*H + K1*K2*K3)
	%fH2A =(C_a0*H^2*K1)/(H^3 + K1*H^2 + K1*K2*H + K1*K2*K3)
	%fH3A =(C_a0*H^3)/(H^3 + K1*H^2 + K1*K2*H + K1*K2*K3)
	%OH =Kw/H

	H3A(i,1) = (C_a0*H(i,1)*H(i,1)^3)/(H(i,1)^3 + K1*H(i,1)^2 + K1*K2*H(i,1) + K1*K2*K3);
	H2A(i,1) = (C_a0*H(i,1)*K1*H(i,1)^2)/(H(i,1)^3 + K1*H(i,1)^2 + K1*K2*H(i,1) + K1*K2*K3);
	HA(i,1) = (C_a0*H(i,1)*K1*K2)/(H(i,1)^3 + K1*H(i,1)^2 + K1*K2*H(i,1) + K1*K2*K3);
	A(i,1) = (C_a0*K1*K2*K3)/(H(i,1)^3 + K1*H(i,1)^2 + K1*K2*H(i,1) + K1*K2*K3);
	OH(i,1) = Kw/H(i);
	FH3A(i,1)=C_a0;
	FNa(i,1)=C_b0;
	if(i>1)
	Beta(i,1)=(FNa(i,1)-FNa(i-1,1))/(PH(i,1)-PH(i-1,1));
	else 
	Beta(i,1)=0;
	end

	fprintf(outfile,'%.12f %.12f %.12f %.12f %.12f %.12f %.12f %.20f %.20f %.12f %.12f %.12f %.12f\n',Acid_vol, x(i), vol_tot, H3A(i,1), H2A(i,1),HA(i,1),A(i,1),H(i,1), OH(i,1), PH(i,1),FH3A(i,1),FNa(i,1),Beta(i,1));

	

end
fclose(outfile);


% Plot it
plot(x,PH,'red')
xlabel('Volume Alkali (ml)')
ylabel('pH')
title('Titration (Volume of alkali vs pH')
