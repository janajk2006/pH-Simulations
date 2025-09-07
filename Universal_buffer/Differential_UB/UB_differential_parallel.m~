% Obtained form the webpage: http://www.chemicalforums.com/index.php?topic=80995.0
%Takes 850 sec on 8 cores, 8 GB memory Dell machine

clear
clc
tic

Rate_const_no=21;
Diffeq_no=15;

%Numerical integration
Custom_RelTol=1e-10;
Custom_AbsTol=1e-18;
Max_time=1e18; %sec


%Cluster property
processor_no=8;
Cluster_property = parcluster;
Cluster_property.NumWorkers = 8;
saveProfile(Cluster_property);

% Equilibrium constants
KW = power(10,-14);
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
C_acid_a = 1; % Molarity
C_acid_b = 1; % Molarity
C_acid_c = 1; % Molarity
C_base = 1; % 1.0

% Initial volume of acid
Acid_vol_a = 10; % ml
Acid_vol_b = 10; % ml
Acid_vol_c = 10; % ml
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
NAOH=size(Titration_resolution,1);
H=size(Titration_resolution,1);
OH=size(Titration_resolution,1);
PH=size(Titration_resolution,1);

H3A_vol=size(Titration_resolution,1);
H3B_vol=size(Titration_resolution,1);
H3C_vol=size(Titration_resolution,1);
NAOH_vol=size(Titration_resolution,1);
Tot_vol=size(Titration_resolution,1);


FH3A=size(Titration_resolution,1);	% Sodium phosphate concentration
FH3B=size(Titration_resolution,1);	% Citric acid concentration
FH3C=size(Titration_resolution,1);	% Boric acid concentration 
FNa=size(Titration_resolution,1);	% NaOH concentration
Beta=size(Titration_resolution,1);


K=size(Rate_const_no);
Y_init=size(Titration_resolution,Diffeq_no);
Y_tol=size(Diffeq_no);

for i = 1:Diffeq_no
Y_tol(i)=Custom_AbsTol;
end

outfile=fopen('Results','w');
fprintf(outfile,'AA_Acid_vol AB_Acid_vol AC_Acid_vol  x(i)  vol_tot  H3A(i,1)  H2A(i,1) HA(i,1) A(i,1) H3B(i,1)  H2B(i,1) HB(i,1) B(i,1) H3C(i,1)  H2C(i,1) HC(i,1) C(i,1) H(i,1)  OH(i,1)  PH(i,1) FH3A(i,1) FH3B(i,1) FH3C(i,1) FNa(i,1) Beta(i,1)\n');



for i = 1:Titration_resolution
	% The total volume
	vol_tot = Acid_vol_a + Acid_vol_b + Acid_vol_c + x(i);

	% Concentrations of acid and base
	C_a0 = (C_acid_a * Acid_vol_a) / vol_tot;
	C_b0 = (C_acid_b * Acid_vol_b) / vol_tot;
	C_c0 = (C_acid_c * Acid_vol_c) / vol_tot;
	C_n0 = (C_base * x(i)) / vol_tot ;

	H3A_vol(i,1)=Acid_vol_a;
	H3B_vol(i,1)=Acid_vol_b;
	H3C_vol(i,1)=Acid_vol_c;
	NAOH_vol(i,1)=x(i);
	Tot_vol(i,1)=vol_tot;


	% Compute the concentration of hydrogen ions

		Y_init(i,1)=C_a0;
		Y_init(i,2)=0;
		Y_init(i,3)=0;
		Y_init(i,4)=0;
		Y_init(i,5)=C_b0;
		Y_init(i,6)=0;
		Y_init(i,7)=0;
		Y_init(i,8)=0;
		Y_init(i,9)=C_c0;
		Y_init(i,10)=0;
		Y_init(i,11)=0;
		Y_init(i,12)=0;
		Y_init(i,13)=power(10,-7);
		Y_init(i,14)=power(10,-7);
		Y_init(i,15)=C_n0;
	
		FH3A(i,1)=C_a0;
		FH3B(i,1)=C_b0;
		FH3C(i,1)=C_c0;
		FNa(i,1)=C_n0;


end

		options = odeset('RelTol',Custom_RelTol,'AbsTol',Y_tol);

		K(1)=K1;		%1.7378*power(10,-5);
		K(2)=power(10,0);
		K(3)=K2;		%6.3*power(10,-1);
		K(4)=power(10,0);	%1.7378*power(10,-5);
		K(5)=K3;		%6.3*power(10,-1);
		K(6)=power(10,0);	%1.7378*power(10,-5);
		K(7)=K4;
		K(8)=power(10,0);
		K(9)=K5;
		K(10)=power(10,0);	%1.7378*power(10,-5);
		K(11)=K6;
		K(12)=power(10,0);		%6.3*power(10,-1);
		K(13)=K7;	%1.7378*power(10,-5);
		K(14)=power(10,0);		%6.3*power(10,-1);
		K(15)=K8;	%1.7378*power(10,-5);
		K(16)=power(10,0);
		K(17)=K9;
		K(18)=power(10,0);
		K(19)=KW;		%1.7378*power(10,-5);
		K(20)=power(10,0);
		K(21)=power(10,2);



matlabpool ('open',processor_no);
parfor i = 1:Titration_resolution


	[T,Y] = ode15s(@UB_ODE,[0,Max_time],Y_init(i,:),options,K);

	%Y = zeros(Titration_resolution,15); 
	S=size(Y);

	% Calculate all other species
	H3A(i,1) = Y(S(1),1);
	H2A(i,1) = Y(S(1),2);
	HA(i,1) = Y(S(1),3);
	A(i,1) = Y(S(1),4);
	H3B(i,1) = Y(S(1),5);
	H2B(i,1) = Y(S(1),6);
	HB(i,1) = Y(S(1),7);
	B(i,1) = Y(S(1),8);
	H3C(i,1) = Y(S(1),9);
	H2C(i,1) = Y(S(1),10);
	HC(i,1) = Y(S(1),11);
	C(i,1) = Y(S(1),12);
	H(i,1) = Y(S(1),13);
	OH(i,1) = Y(S(1),14);
	NAOH(i,1) = Y(S(1),15);

	% Calculate the pH value
	PH(i,1) = -log10(H(i,1));	

i
end


matlabpool close;


	H3A_vol(i,1)=Acid_vol_a;
	H3B_vol(i,1)=Acid_vol_b;
	H3C_vol(i,1)=Acid_vol_c;
	NAOH_vol(i,1)=x(i);
	Tot_vol(i,1)=vol_tot;

% Print the results
for i = 1:Titration_resolution

	if(i>1)
	Beta(i,1)=(FNa(i,1)-FNa(i-1,1))/(PH(i,1)-PH(i-1,1));
	else 
	Beta(i,1)=0;
	end


fprintf(outfile,'%.12f %.12f %.12f %.12f %.12f %.12f %.12f %.12f %.12f %.12f %.12f %.12f %.12f %.12f %.12f %.12f %.12f %.20f %.20f %.12f %.12f %.12f %.12f %.12f %.12f\n',H3A_vol(i),H3B_vol(i),H3C_vol(i), x(i), Tot_vol(i), H3A(i,1), H2A(i,1),HA(i,1),A(i,1),H3B(i,1), H2B(i,1),HB(i,1),B(i,1),H3C(i,1), H2C(i,1),HC(i,1),C(i,1),H(i,1), OH(i,1), PH(i,1),FH3A(i,1),FH3B(i,1),FH3C(i,1),FNa(i,1),Beta(i,1));

end

fclose(outfile);


% Plot it
figure();
plot(x,PH,'red')
xlabel('Volume Alkali (ml)')
ylabel('pH')
title('Titration (Volume of alkali vs pH')
toc
