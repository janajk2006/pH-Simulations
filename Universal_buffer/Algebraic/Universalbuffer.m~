% Obtained form the webpage: http://www.chemicalforums.com/index.php?topic=80995.0
clear;
clc;

Max_species=9;
Max_roots=Max_species+2;
% Equilibrium constants
Kw = power(10,-14);
K1 = power(10,-2.12);		%power(10,-3.08);
K2 = power(10,-7.21); 		%power(10,-4.74);
K3 = power(10,-12.32); 		%power(10,-5.4);
K4 = power(10,-3.08);		%power(10,-3.08);
K5 = power(10,-4.74); 		%power(10,-4.74);
K6 = power(10,-5.4); 		%power(10,-5.4);
K7 = power(10,-9.23);		%power(10,-3.08);
K8 = power(10,-12.4); 		%power(10,-4.74);
K9 = power(10,-13.3); 		%power(10,-5.4);

% Initial concentration of acid
AA_C_acid = 1; % mol/dm^3
AB_C_acid = 1; % mol/dm^3
AC_C_acid = 1; % mol/dm^3
C_base = 1; % mol/dm^3

% Initial volume of acid
AA_Acid_vol = 10; % ml
AB_Acid_vol = 10; % ml
AC_Acid_vol = 10; % ml
Base_vol = 150; % ml



% Graph settings
Titration_resolution=500; % number of titration points
dx = Base_vol/(Titration_resolution-1);
x = 0.0:dx:Base_vol;


% Array for all the species
H3A=size(Titration_resolution,1);
H2A=size(Titration_resolution,1);
HA=size(Titration_resolution,1);
A=size(Titration_resolution,1);
H3B=size(Titration_resolution,1);
H2B=size(Titration_resolution,1);
HB=size(Titration_resolution,1);
B=size(Titration_resolution,1);
H3C=size(Titration_resolution,1);
H2C=size(Titration_resolution,1);
HC=size(Titration_resolution,1);
C=size(Titration_resolution,1);
H=size(Titration_resolution,1);
PH=size(Titration_resolution,1);
OH=size(Titration_resolution,1);
FNa=size(Titration_resolution,1);
FH3A=size(Titration_resolution,1);
FH3B=size(Titration_resolution,1);
FH3C=size(Titration_resolution,1);
Beta=size(Titration_resolution,1);
outfile=fopen('Results','w');

for j = 1:Max_roots
fprintf(outfile,'Roots(%d) ',j);
end
fprintf(outfile,'AA_Acid_vol AB_Acid_vol AC_Acid_vol  x(i)  vol_tot  H3A(i,1)  H2A(i,1) HA(i,1) A(i,1) H3B(i,1)  H2B(i,1) HB(i,1) B(i,1) H3C(i,1)  H2C(i,1) HC(i,1) C(i,1) H(i,1)  OH(i,1)  PH(i,1) FH3A(i,1) FH3B(i,1) FH3C(i,1) FNa(i,1) Beta(i,1)\n');


for i = 1:Titration_resolution

	% The total volume
	acid_tot=AA_Acid_vol + AB_Acid_vol + AC_Acid_vol;
	vol_tot = acid_tot + x(i);

	% Concentrations of acid and base
	AA_C_a0 = (AA_C_acid * AA_Acid_vol) / vol_tot;
	AB_C_a0 = (AB_C_acid * AB_Acid_vol) / vol_tot;
	AC_C_a0 = (AC_C_acid * AC_Acid_vol) / vol_tot;
	C_b0 = (C_base * x(i)) / vol_tot ;


	% Compute the concentration of hydrogen ions
	% Solve the polynomial equation
	%


	R = roots([1, (C_b0 + K1 + K4 + K7), (C_b0*K1 - AA_C_a0*K1 - AB_C_a0*K4 - AC_C_a0*K7 - Kw + C_b0*K4 + C_b0*K7 + K1*K2 + K1*K4 + K1*K7 + K4*K5 + K4*K7 + K7*K8) , (C_b0*K1*K2 - K4*Kw - K7*Kw - 2*AA_C_a0*K1*K2 - AA_C_a0*K1*K4 - AA_C_a0*K1*K7 - AB_C_a0*K1*K4 - 2*AB_C_a0*K4*K5 - AB_C_a0*K4*K7 - AC_C_a0*K1*K7 - AC_C_a0*K4*K7 - 2*AC_C_a0*K7*K8 - K1*Kw + C_b0*K1*K4 + C_b0*K1*K7 + C_b0*K4*K5 + C_b0*K4*K7 + C_b0*K7*K8 + K1*K2*K3 + K1*K2*K4 + K1*K2*K7 + K1*K4*K5 + K1*K4*K7 + K4*K5*K6 + K1*K7*K8 + K4*K5*K7 + K4*K7*K8 + K7*K8*K9) , (C_b0*K1*K2*K3 - K1*K4*Kw - K1*K7*Kw - K4*K5*Kw - K4*K7*Kw - K7*K8*Kw - 3*AA_C_a0*K1*K2*K3 - 2*AA_C_a0*K1*K2*K4 - 2*AA_C_a0*K1*K2*K7 - AA_C_a0*K1*K4*K5 - AA_C_a0*K1*K4*K7 - AA_C_a0*K1*K7*K8 - AB_C_a0*K1*K2*K4 - 2*AB_C_a0*K1*K4*K5 - AB_C_a0*K1*K4*K7 - 3*AB_C_a0*K4*K5*K6 - 2*AB_C_a0*K4*K5*K7 - AB_C_a0*K4*K7*K8 - AC_C_a0*K1*K2*K7 - AC_C_a0*K1*K4*K7 - 2*AC_C_a0*K1*K7*K8 - AC_C_a0*K4*K5*K7 - 2*AC_C_a0*K4*K7*K8 - 3*AC_C_a0*K7*K8*K9 - K1*K2*Kw + C_b0*K1*K2*K4 + C_b0*K1*K2*K7 + C_b0*K1*K4*K5 + C_b0*K1*K4*K7 + C_b0*K4*K5*K6 + C_b0*K1*K7*K8 + C_b0*K4*K5*K7 + C_b0*K4*K7*K8 + C_b0*K7*K8*K9 + K1*K2*K3*K4 + K1*K2*K4*K5 + K1*K2*K3*K7 + K1*K2*K4*K7 + K1*K4*K5*K6 + K1*K4*K5*K7 + K1*K2*K7*K8 + K1*K4*K7*K8 + K4*K5*K6*K7 + K4*K5*K7*K8 + K1*K7*K8*K9 + K4*K7*K8*K9) , (C_b0*K1*K2*K3*K4 - K1*K2*K4*Kw - K1*K2*K7*Kw - K1*K4*K5*Kw - K1*K4*K7*Kw - K4*K5*K6*Kw - K1*K7*K8*Kw - K4*K5*K7*Kw - K4*K7*K8*Kw - K7*K8*K9*Kw - 3*AA_C_a0*K1*K2*K3*K4 - 2*AA_C_a0*K1*K2*K4*K5 - 3*AA_C_a0*K1*K2*K3*K7 - 2*AA_C_a0*K1*K2*K4*K7 - AA_C_a0*K1*K4*K5*K6 - AA_C_a0*K1*K4*K5*K7 - 2*AA_C_a0*K1*K2*K7*K8 - AA_C_a0*K1*K4*K7*K8 - AA_C_a0*K1*K7*K8*K9 - AB_C_a0*K1*K2*K3*K4 - 2*AB_C_a0*K1*K2*K4*K5 - AB_C_a0*K1*K2*K4*K7 - 3*AB_C_a0*K1*K4*K5*K6 - 2*AB_C_a0*K1*K4*K5*K7 - AB_C_a0*K1*K4*K7*K8 - 3*AB_C_a0*K4*K5*K6*K7 - 2*AB_C_a0*K4*K5*K7*K8 - AB_C_a0*K4*K7*K8*K9 - AC_C_a0*K1*K2*K3*K7 - AC_C_a0*K1*K2*K4*K7 - AC_C_a0*K1*K4*K5*K7 - 2*AC_C_a0*K1*K2*K7*K8 - 2*AC_C_a0*K1*K4*K7*K8 - AC_C_a0*K4*K5*K6*K7 - 2*AC_C_a0*K4*K5*K7*K8 - 3*AC_C_a0*K1*K7*K8*K9 - 3*AC_C_a0*K4*K7*K8*K9 - K1*K2*K3*Kw + C_b0*K1*K2*K4*K5 + C_b0*K1*K2*K3*K7 + C_b0*K1*K2*K4*K7 + C_b0*K1*K4*K5*K6 + C_b0*K1*K4*K5*K7 + C_b0*K1*K2*K7*K8 + C_b0*K1*K4*K7*K8 + C_b0*K4*K5*K6*K7 + C_b0*K4*K5*K7*K8 + C_b0*K1*K7*K8*K9 + C_b0*K4*K7*K8*K9 + K1*K2*K3*K4*K5 + K1*K2*K3*K4*K7 + K1*K2*K4*K5*K6 + K1*K2*K4*K5*K7 + K1*K2*K3*K7*K8 + K1*K2*K4*K7*K8 + K1*K4*K5*K6*K7 + K1*K4*K5*K7*K8 + K1*K2*K7*K8*K9 + K1*K4*K7*K8*K9 + K4*K5*K6*K7*K8 + K4*K5*K7*K8*K9) , (C_b0*K1*K2*K3*K4*K5 - K1*K2*K4*K5*Kw - K1*K2*K3*K7*Kw - K1*K2*K4*K7*Kw - K1*K4*K5*K6*Kw - K1*K4*K5*K7*Kw - K1*K2*K7*K8*Kw - K1*K4*K7*K8*Kw - K4*K5*K6*K7*Kw - K4*K5*K7*K8*Kw - K1*K7*K8*K9*Kw - K4*K7*K8*K9*Kw - 3*AA_C_a0*K1*K2*K3*K4*K5 - 3*AA_C_a0*K1*K2*K3*K4*K7 - 2*AA_C_a0*K1*K2*K4*K5*K6 - 2*AA_C_a0*K1*K2*K4*K5*K7 - 3*AA_C_a0*K1*K2*K3*K7*K8 - 2*AA_C_a0*K1*K2*K4*K7*K8 - AA_C_a0*K1*K4*K5*K6*K7 - AA_C_a0*K1*K4*K5*K7*K8 - 2*AA_C_a0*K1*K2*K7*K8*K9 - AA_C_a0*K1*K4*K7*K8*K9 - 2*AB_C_a0*K1*K2*K3*K4*K5 - AB_C_a0*K1*K2*K3*K4*K7 - 3*AB_C_a0*K1*K2*K4*K5*K6 - 2*AB_C_a0*K1*K2*K4*K5*K7 - AB_C_a0*K1*K2*K4*K7*K8 - 3*AB_C_a0*K1*K4*K5*K6*K7 - 2*AB_C_a0*K1*K4*K5*K7*K8 - AB_C_a0*K1*K4*K7*K8*K9 - 3*AB_C_a0*K4*K5*K6*K7*K8 - 2*AB_C_a0*K4*K5*K7*K8*K9 - AC_C_a0*K1*K2*K3*K4*K7 - AC_C_a0*K1*K2*K4*K5*K7 - 2*AC_C_a0*K1*K2*K3*K7*K8 - 2*AC_C_a0*K1*K2*K4*K7*K8 - AC_C_a0*K1*K4*K5*K6*K7 - 2*AC_C_a0*K1*K4*K5*K7*K8 - 3*AC_C_a0*K1*K2*K7*K8*K9 - 3*AC_C_a0*K1*K4*K7*K8*K9 - 2*AC_C_a0*K4*K5*K6*K7*K8 - 3*AC_C_a0*K4*K5*K7*K8*K9 - K1*K2*K3*K4*Kw + C_b0*K1*K2*K3*K4*K7 + C_b0*K1*K2*K4*K5*K6 + C_b0*K1*K2*K4*K5*K7 + C_b0*K1*K2*K3*K7*K8 + C_b0*K1*K2*K4*K7*K8 + C_b0*K1*K4*K5*K6*K7 + C_b0*K1*K4*K5*K7*K8 + C_b0*K1*K2*K7*K8*K9 + C_b0*K1*K4*K7*K8*K9 + C_b0*K4*K5*K6*K7*K8 + C_b0*K4*K5*K7*K8*K9 + K1*K2*K3*K4*K5*K6 + K1*K2*K3*K4*K5*K7 + K1*K2*K3*K4*K7*K8 + K1*K2*K4*K5*K6*K7 + K1*K2*K4*K5*K7*K8 + K1*K2*K3*K7*K8*K9 + K1*K2*K4*K7*K8*K9 + K1*K4*K5*K6*K7*K8 + K1*K4*K5*K7*K8*K9 + K4*K5*K6*K7*K8*K9), (C_b0*K1*K2*K3*K4*K5*K6 - K1*K2*K3*K4*K7*Kw - K1*K2*K4*K5*K6*Kw - K1*K2*K4*K5*K7*Kw - K1*K2*K3*K7*K8*Kw - K1*K2*K4*K7*K8*Kw - K1*K4*K5*K6*K7*Kw - K1*K4*K5*K7*K8*Kw - K1*K2*K7*K8*K9*Kw - K1*K4*K7*K8*K9*Kw - K4*K5*K6*K7*K8*Kw - K4*K5*K7*K8*K9*Kw - 3*AA_C_a0*K1*K2*K3*K4*K5*K6 - 3*AA_C_a0*K1*K2*K3*K4*K5*K7 - 3*AA_C_a0*K1*K2*K3*K4*K7*K8 - 2*AA_C_a0*K1*K2*K4*K5*K6*K7 - 2*AA_C_a0*K1*K2*K4*K5*K7*K8 - 3*AA_C_a0*K1*K2*K3*K7*K8*K9 - 2*AA_C_a0*K1*K2*K4*K7*K8*K9 - AA_C_a0*K1*K4*K5*K6*K7*K8 - AA_C_a0*K1*K4*K5*K7*K8*K9 - 3*AB_C_a0*K1*K2*K3*K4*K5*K6 - 2*AB_C_a0*K1*K2*K3*K4*K5*K7 - AB_C_a0*K1*K2*K3*K4*K7*K8 - 3*AB_C_a0*K1*K2*K4*K5*K6*K7 - 2*AB_C_a0*K1*K2*K4*K5*K7*K8 - AB_C_a0*K1*K2*K4*K7*K8*K9 - 3*AB_C_a0*K1*K4*K5*K6*K7*K8 - 2*AB_C_a0*K1*K4*K5*K7*K8*K9 - 3*AB_C_a0*K4*K5*K6*K7*K8*K9 - AC_C_a0*K1*K2*K3*K4*K5*K7 - 2*AC_C_a0*K1*K2*K3*K4*K7*K8 - AC_C_a0*K1*K2*K4*K5*K6*K7 - 2*AC_C_a0*K1*K2*K4*K5*K7*K8 - 3*AC_C_a0*K1*K2*K3*K7*K8*K9 - 3*AC_C_a0*K1*K2*K4*K7*K8*K9 - 2*AC_C_a0*K1*K4*K5*K6*K7*K8 - 3*AC_C_a0*K1*K4*K5*K7*K8*K9 - 3*AC_C_a0*K4*K5*K6*K7*K8*K9 - K1*K2*K3*K4*K5*Kw + C_b0*K1*K2*K3*K4*K5*K7 + C_b0*K1*K2*K3*K4*K7*K8 + C_b0*K1*K2*K4*K5*K6*K7 + C_b0*K1*K2*K4*K5*K7*K8 + C_b0*K1*K2*K3*K7*K8*K9 + C_b0*K1*K2*K4*K7*K8*K9 + C_b0*K1*K4*K5*K6*K7*K8 + C_b0*K1*K4*K5*K7*K8*K9 + C_b0*K4*K5*K6*K7*K8*K9 + K1*K2*K3*K4*K5*K6*K7 + K1*K2*K3*K4*K5*K7*K8 + K1*K2*K4*K5*K6*K7*K8 + K1*K2*K3*K4*K7*K8*K9 + K1*K2*K4*K5*K7*K8*K9 + K1*K4*K5*K6*K7*K8*K9) , (C_b0*K1*K2*K3*K4*K5*K6*K7 - K1*K2*K3*K4*K5*K7*Kw - K1*K2*K3*K4*K7*K8*Kw - K1*K2*K4*K5*K6*K7*Kw - K1*K2*K4*K5*K7*K8*Kw - K1*K2*K3*K7*K8*K9*Kw - K1*K2*K4*K7*K8*K9*Kw - K1*K4*K5*K6*K7*K8*Kw - K1*K4*K5*K7*K8*K9*Kw - K4*K5*K6*K7*K8*K9*Kw - 3*AA_C_a0*K1*K2*K3*K4*K5*K6*K7 - 3*AA_C_a0*K1*K2*K3*K4*K5*K7*K8 - 2*AA_C_a0*K1*K2*K4*K5*K6*K7*K8 - 3*AA_C_a0*K1*K2*K3*K4*K7*K8*K9 - 2*AA_C_a0*K1*K2*K4*K5*K7*K8*K9 - AA_C_a0*K1*K4*K5*K6*K7*K8*K9 - 3*AB_C_a0*K1*K2*K3*K4*K5*K6*K7 - 2*AB_C_a0*K1*K2*K3*K4*K5*K7*K8 - 3*AB_C_a0*K1*K2*K4*K5*K6*K7*K8 - AB_C_a0*K1*K2*K3*K4*K7*K8*K9 - 2*AB_C_a0*K1*K2*K4*K5*K7*K8*K9 - 3*AB_C_a0*K1*K4*K5*K6*K7*K8*K9 - AC_C_a0*K1*K2*K3*K4*K5*K6*K7 - 2*AC_C_a0*K1*K2*K3*K4*K5*K7*K8 - 2*AC_C_a0*K1*K2*K4*K5*K6*K7*K8 - 3*AC_C_a0*K1*K2*K3*K4*K7*K8*K9 - 3*AC_C_a0*K1*K2*K4*K5*K7*K8*K9 - 3*AC_C_a0*K1*K4*K5*K6*K7*K8*K9 - K1*K2*K3*K4*K5*K6*Kw + C_b0*K1*K2*K3*K4*K5*K7*K8 + C_b0*K1*K2*K4*K5*K6*K7*K8 + C_b0*K1*K2*K3*K4*K7*K8*K9 + C_b0*K1*K2*K4*K5*K7*K8*K9 + C_b0*K1*K4*K5*K6*K7*K8*K9 + K1*K2*K3*K4*K5*K6*K7*K8 + K1*K2*K3*K4*K5*K7*K8*K9 + K1*K2*K4*K5*K6*K7*K8*K9), (K1*K2*K3*K4*K5*K6*K7*K8*K9 - K1*K2*K3*K4*K5*K7*K8*Kw - K1*K2*K4*K5*K6*K7*K8*Kw - K1*K2*K3*K4*K7*K8*K9*Kw - K1*K2*K4*K5*K7*K8*K9*Kw - K1*K4*K5*K6*K7*K8*K9*Kw - K1*K2*K3*K4*K5*K6*K7*Kw - 3*AA_C_a0*K1*K2*K3*K4*K5*K6*K7*K8 - 3*AA_C_a0*K1*K2*K3*K4*K5*K7*K8*K9 - 2*AA_C_a0*K1*K2*K4*K5*K6*K7*K8*K9 - 3*AB_C_a0*K1*K2*K3*K4*K5*K6*K7*K8 - 2*AB_C_a0*K1*K2*K3*K4*K5*K7*K8*K9 - 3*AB_C_a0*K1*K2*K4*K5*K6*K7*K8*K9 - 2*AC_C_a0*K1*K2*K3*K4*K5*K6*K7*K8 - 3*AC_C_a0*K1*K2*K3*K4*K5*K7*K8*K9 - 3*AC_C_a0*K1*K2*K4*K5*K6*K7*K8*K9 + C_b0*K1*K2*K3*K4*K5*K6*K7*K8 + C_b0*K1*K2*K3*K4*K5*K7*K8*K9 + C_b0*K1*K2*K4*K5*K6*K7*K8*K9), (C_b0*K1*K2*K3*K4*K5*K6*K7*K8*K9 - K1*K2*K3*K4*K5*K7*K8*K9*Kw - K1*K2*K4*K5*K6*K7*K8*K9*Kw - 3*AA_C_a0*K1*K2*K3*K4*K5*K6*K7*K8*K9 - 3*AB_C_a0*K1*K2*K3*K4*K5*K6*K7*K8*K9 - 3*AC_C_a0*K1*K2*K3*K4*K5*K6*K7*K8*K9 - K1*K2*K3*K4*K5*K6*K7*K8*Kw),- K1*K2*K3*K4*K5*K6*K7*K8*K9*Kw]);

	% Remove all chemically irrelevant solutions	
	check_index=1;
	Temp_H=0;
	for j = 1:Max_roots
		fprintf(outfile,'%.20f ',R(j,1));
		if R(j,1)==real(R(j,1)) && R(j,1)>0 && R(j,1)<(AA_C_a0 + AB_C_a0 + AC_C_a0)
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

	%fA =(AA_C_a0*K1*K2*K3)/(H^3 + K1*H^2 + K1*K2*H + K1*K2*K3) 
	%fHA =(AA_C_a0*H*K1*K2)/(H^3 + K1*H^2 + K1*K2*H + K1*K2*K3) 
	%fH2A =(AA_C_a0*H^2*K1)/(H^3 + K1*H^2 + K1*K2*H + K1*K2*K3)
	%fH3A =(AA_C_a0*H^3)/(H^3 + K1*H^2 + K1*K2*H + K1*K2*K3)
	 
	%fB =(AB_C_a0*K4*K5*K6)/(H^3 + K4*H^2 + K4*K5*H + K4*K5*K6) 
	%fHB =(AB_C_a0*H*K4*K5)/(H^3 + K4*H^2 + K4*K5*H + K4*K5*K6) 
	%fH2B =(AB_C_a0*H^2*K4)/(H^3 + K4*H^2 + K4*K5*H + K4*K5*K6) 
	%fH3B =(AB_C_a0*H^3)/(H^3 + K4*H^2 + K4*K5*H + K4*K5*K6)
	 
	%fC =(AC_C_a0*K7*K8*K9)/(H^3 + K7*H^2 + K7*K8*H + K7*K8*K9) 
	%fHC =(AC_C_a0*H*K7*K8)/(H^3 + K7*H^2 + K7*K8*H + K7*K8*K9)
	%fH2C =(AC_C_a0*H^2*K7)/(H^3 + K7*H^2 + K7*K8*H + K7*K8*K9) 
	%fH3C =(AC_C_a0*H^3)/(H^3 + K7*H^2 + K7*K8*H + K7*K8*K9)
	 
	%W =Kw/H

	H3A(i,1) = (AA_C_a0*H(i,1)*H(i,1)^3)/(H(i,1)^3 + K1*H(i,1)^2 + K1*K2*H(i,1) + K1*K2*K3);
	H2A(i,1) = (AA_C_a0*H(i,1)*K1*H(i,1)^2)/(H(i,1)^3 + K1*H(i,1)^2 + K1*K2*H(i,1) + K1*K2*K3);
	HA(i,1) = (AA_C_a0*H(i,1)*K1*K2)/(H(i,1)^3 + K1*H(i,1)^2 + K1*K2*H(i,1) + K1*K2*K3);
	A(i,1) = (AA_C_a0*K1*K2*K3)/(H(i,1)^3 + K1*H(i,1)^2 + K1*K2*H(i,1) + K1*K2*K3);

	H3B(i,1) =(AB_C_a0*K4*K5*K6)/(H(i,1)^3 + K4*H(i,1)^2 + K4*K5*H(i,1) + K4*K5*K6);
	H2B(i,1) =(AB_C_a0*H(i,1)*K4*K5)/(H(i,1)^3 + K4*H(i,1)^2 + K4*K5*H(i,1) + K4*K5*K6);
	HB(i,1) =(AB_C_a0*H(i,1)^2*K4)/(H(i,1)^3 + K4*H(i,1)^2 + K4*K5*H(i,1) + K4*K5*K6);
	B(i,1) =(AB_C_a0*H(i,1)^3)/(H(i,1)^3 + K4*H(i,1)^2 + K4*K5*H(i,1) + K4*K5*K6);

	H3C(i,1) =(AC_C_a0*K7*K8*K9)/(H(i,1)^3 + K7*H(i,1)^2 + K7*K8*H(i,1) + K7*K8*K9); 
	H2C(i,1) =(AC_C_a0*H(i,1)*K7*K8)/(H(i,1)^3 + K7*H(i,1)^2 + K7*K8*H(i,1) + K7*K8*K9);
	HC(i,1) =(AC_C_a0*H(i,1)^2*K7)/(H(i,1)^3 + K7*H(i,1)^2 + K7*K8*H(i,1) + K7*K8*K9);
	C(i,1) =(AC_C_a0*H(i,1)^3)/(H(i,1)^3 + K7*H(i,1)^2 + K7*K8*H(i,1) + K7*K8*K9);

	OH(i,1) = Kw/H(i);
	FH3A(i,1)=AA_C_a0;
	FH3B(i,1)=AB_C_a0;
	FH3C(i,1)=AC_C_a0;
	FNa(i,1)=C_b0;
	if(i>1)
	Beta(i,1)=(FNa(i,1)-FNa(i-1,1))/(PH(i,1)-PH(i-1,1));
	else 
	Beta(i,1)=0;
	end

	fprintf(outfile,'%.12f %.12f %.12f %.12f %.12f %.12f %.12f %.12f %.12f %.12f %.12f %.12f %.12f %.12f %.12f %.12f %.12f %.20f %.20f %.12f %.12f %.12f %.12f %.12f %.12f\n',AA_Acid_vol,AB_Acid_vol,AC_Acid_vol, x(i), vol_tot, H3A(i,1), H2A(i,1),HA(i,1),A(i,1),H3B(i,1), H2B(i,1),HB(i,1),B(i,1),H3C(i,1), H2C(i,1),HC(i,1),C(i,1),H(i,1), OH(i,1), PH(i,1),FH3A(i,1),FH3B(i,1),FH3C(i,1),FNa(i,1),Beta(i,1));

	

end
fclose(outfile);


% Plot it
figure();
plot(x,PH,'red')
xlabel('Volume Alkali (ml)')
ylabel('pH')
title('Titration (Volume of alkali vs pH)')
