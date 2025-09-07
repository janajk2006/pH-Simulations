% Obtained form the webpage: http://www.chemicalforums.com/index.php?topic=80995.0
% Takes  sec to complete on 8 cores, 8 GB memory Dell machine

clear
clc
tic


Rate_const_no=7;
Diffeq_no=6;

%Numerical integration
Custom_RelTol=1e-10;
Custom_AbsTol=1e-18;
Max_time=1e18; %sec


% Equilibrium constants
KW = power(10,-14);
K1 = power(10,-4);		%power(10,-2);
K2 = power(10,-9); 		%power(10,-5);

% Initial concentration of acid
C_acid = 1; % mol/dm^3
C_base = 1; % mol/dm^3

% Initial volume of acid
Acid_vol = 10; % ml
Base_vol = 35; % ml



% Graph settings
Titration_resolution=500; % number of titration points
dx = Base_vol/(Titration_resolution-1);
x = 0.0:dx:Base_vol;


% Array for all the species
H2A=size(Titration_resolution,1);
HA=size(Titration_resolution,1);
A=size(Titration_resolution,1);
H=size(Titration_resolution,1);
PH=size(Titration_resolution,1);
NAOH=size(Titration_resolution,1);
OH=size(Titration_resolution,1);
FNa=size(Titration_resolution,1);
FH2A=size(Titration_resolution,1);


% Volume and Beta
H2A_vol=size(Titration_resolution,1);
NaOH_vol=size(Titration_resolution,1);
Tot_vol=size(Titration_resolution,1);
Beta=size(Titration_resolution,1);

K=size(Rate_const_no);
Y_init=size(Diffeq_no);

for i = 1:Diffeq_no
Y_tol(i)=Custom_AbsTol;
end


outfile=fopen('Results','w');
fprintf(outfile,'Acid_vol Base_vol Vol_tot H2A(i,1) HA(i,1) A(i,1) H(i,1) OH(i,1) PH(i,1) FH2A(i,1) FNa(i,1) Beta(i,1)\n');


for i = 1:Titration_resolution

	% The total volume
	vol_tot = Acid_vol + x(i);

	% Concentrations of acid and base
	C_a0 = (C_acid * Acid_vol) / vol_tot;
	C_b0 = (C_base * x(i)) / vol_tot ;



	H2A_vol(i,1)=Acid_vol;
	NaOH_vol(i,1)=x(i);
	Tot_vol(i,1)=vol_tot;

	% Compute the concentration of hydrogen ions

		Y_init(1)=C_a0;
		Y_init(2)=0;
		Y_init(3)=0;
		Y_init(4)=C_b0;
		Y_init(5)=power(10,-7);
		Y_init(6)=power(10,-7);


		options = odeset('RelTol',Custom_RelTol,'AbsTol',Y_tol);



		K(1)=K1;		%1.7378*power(10,-5);
		K(2)=power(10,0);
		K(3)=K2;		%6.3*power(10,-1);
		K(4)=power(10,0);	%1.7378*power(10,-5);
		K(5)=power(10,2);
		K(6)=KW;
		K(7)=power(10,0);

		[T,Y] = ode15s(@Diprotic_ODE,[0,Max_time],Y_init,options,K);


	S=size(Y);




	% Calculate all other species
	H2A(i,1) = Y(S(1),1);
	HA(i,1) = Y(S(1),2);
	A(i,1) = Y(S(1),3);
	NAOH(i,1) = Y(S(1),4);
	H(i,1) = Y(S(1),5);
	OH(i,1) = Y(S(1),6);
	FH2A(i,1)=C_a0;
	FNa(i,1)=C_b0;

	% Calculate the pH value
	PH(i,1) = -log10(H(i,1));	



	if(i>1)
	Beta(i,1)=(FNa(i,1)-FNa(i-1,1))/(PH(i,1)-PH(i-1,1));
	else 
	Beta(i,1)=0;
	end

	fprintf(outfile,'%.12f %.12f %.12f %.12f %.12f %.12f %.20f %.20f %.12f %.12f %.12f %.12f\n',H2A_vol(i,1), NaOH_vol(i,1), Tot_vol(i,1), H2A(i,1),HA(i,1),A(i,1),H(i,1), OH(i,1), PH(i,1),FH2A(i,1),FNa(i,1),Beta(i,1));

i
end

fclose(outfile);
% Plot it
figure();
plot(x,PH,'red')
xlabel('Volume Alkali (ml)')
ylabel('pH')
title('Titration (Volume of alkali vs pH')
toc
