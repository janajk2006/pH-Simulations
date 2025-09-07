

%%%%%%%%%%%%%%%%%%%
% Monoprotic acid %
%%%%%%%%%%%%%%%%%%%

% With Derivation
%----------------

%{	K1
HA <=========> H + A
	Kw'
H20 <=========> H + OH

K1=(H*A)/HA
C_a0=A+HA
Kw=(H*OH)
Na=C_b0
%}

clc;
clear;
syms C_a0 C_b0 Na H OH HA A K1 Kw 
[fA,fHA]=solve(H*A/HA == K1, A+HA == C_a0,A,HA)
W=solve(H*OH==Kw,OH)

%Charge neutrality
%-----------------
X=simplifyFraction(H+C_b0-W-fA,'Expand', true)
[Y_1,Y_2]=numden(X)
Z=collect(Y_1,H)


% Explicitly
%-----------

clc;
clear;
syms C_a0 C_b0 H K1 Kw
X=simplifyFraction(H+C_b0-Kw/H-(K1*C_a0)/(K1+H),'Expand', true)
[Y_1,Y_2]=numden(X)
Z=collect(Y_1,H)

%%%%%%%%%%%%%%%%%
% Diprotic acid %
%%%%%%%%%%%%%%%%%

% With Derivation
%----------------

%{	K1
H2A <=========> HA + H
	K2
HA <=========> A + H
	Kw'
H20 <=========> H + OH

K1=(HA*H)/H2A
K2=(A*H)/HA
C_a0=A+HA+H2A
Kw=(H*OH)
Na=C_b0
%}

clc;
clear;
syms C_a0 C_b0 Na H OH Kw H2A HA A K1 K2
[fA,fHA,fH2A]=solve((HA*H)/H2A == K1, (A*H)/HA == K2, A+HA+H2A == C_a0,A,HA,H2A)
W=solve(H*OH==Kw,OH)

%Charge neutrality
%-----------------
X=simplifyFraction(H+C_b0-W-2*fA-fHA,'Expand', true)
[Y_1,Y_2]=numden(X)
Z=collect(Y_1,H)


% Explicitly
%-----------
clc;
clear;
syms C_a0 C_b0 H K1 K2 Kw OH H2A HA A
X=simplifyFraction(H+C_b0-Kw/H-2*(K1*K2*C_a0)/(H*H+K1*H+K1*K2)-(K1*H*C_a0)/(H*H+K1*H+K1*K2),'Expand', true)
[Y_1,Y_2]=numden(X)
Z=collect(Y_1,H)



%%%%%%%%%%%%%%%%%%
% Triprotic acid %
%%%%%%%%%%%%%%%%%%

% With Derivation
%----------------

%{	K1
H3A <=========> H2A + H
	K2
H2A <=========> HA + H
	K2
HA <=========> A + H

	Kw'
H20 <=========> H + OH

K1=(H2A*H)/H3A
K2=(HA*H)/H2A
K3=(A*H)/HA
C_a0=A+HA+H2A+H3A
Kw=(H*OH)
Na=C_b0
%}

clc;
clear;
syms C_a0 C_b0 Na H OH Kw H3A H2A HA A K1 K2 K3
[fA,fHA,fH2A, fH3A]=solve((H2A*H)/H3A == K1, (HA*H)/H2A == K2, (A*H)/HA == K3, A+HA+H2A+H3A == C_a0,A,HA,H2A,H3A)
W=solve(H*OH==Kw,OH)

%Charge neutrality
%-----------------
X=simplifyFraction(H+C_b0-W-3*fA-2*fHA-fH2A,'Expand', true)
[Y_1,Y_2]=numden(X)
Z=collect(Y_1,H)




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Universal buffer with Triprotic acid 	   %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% With Derivation
%----------------

%{	K1
H3A <=========> H2A + H
	K2
H2A <=========> HA + H
	K3
HA <=========> A + H


	K4
H3B <=========> H2B + H
	K5
H2B <=========> HB + H
	K6
HB <=========> B + H

	K7
H3C <=========> H2C + H
	K8
H2C <=========> HC + H
	K9
HC <=========> C + H


	Kw'
H20 <=========> H + OH

K1=(H2A*H)/H3A
K2=(HA*H)/H2A
K3=(A*H)/HA

K4=(H2B*H)/H3B
K5=(HB*H)/H2B
K6=(B*H)/HB

K7=(H2C*H)/H3C
K8=(HC*H)/H2C
K9=(C*H)/HC

AA_C_a0=A+HA+H2A+H3A
AB_C_a0=B+HB+H2B+H3B
AC_C_a0=C+HC+H2C+H3C
Kw=(H*OH)
Na=C_b0
%}

clc;
clear;
syms AA_C_a0 AB_C_a0 AC_C_a0 C_b0 Na H OH Kw H3A H2A HA A  H3B H2B HB B  H3C H2C HC C K1 K2 K3 K4 K5 K6 K7 K8 K9
[fA,fHA,fH2A,fH3A]=solve((H2A*H)/H3A == K1, (HA*H)/H2A == K2, (A*H)/HA == K3, A+HA+H2A+H3A == AA_C_a0,A,HA,H2A,H3A)

[fB,fHB,fH2B,fH3B]=solve((H2B*H)/H3B == K4, (HB*H)/H2B == K5, (B*H)/HB == K6, B+HB+H2B+H3B == AB_C_a0,B,HB,H2B,H3B)

[fC,fHC,fH2C,fH3C]=solve((H2C*H)/H3C == K7, (HC*H)/H2C == K8, (C*H)/HC == K9, C+HC+H2C+H3C == AC_C_a0,C,HC,H2C,H3C)
W=solve(H*OH==Kw,OH)

%Charge neutrality
%-----------------
X=simplifyFraction(H+C_b0-W-3*fA-2*fHA-fH2A-3*fB-2*fHB-fH2B-3*fC-2*fHC-fH2C,'Expand', true)
[Y_1,Y_2]=numden(X)
Z=collect(Y_1,H)


