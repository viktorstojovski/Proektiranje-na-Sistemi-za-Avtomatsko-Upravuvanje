%% zadaca 3.1
clc
clear
s=tf('s')
P=1/(s*(s+4)^2)
R1=(160*(10*s+1))/(200*s+1)
R2=24*(s+1)/(s+4)
R3=(s+4)/(s+1)
G01=R1*P
G02=R2*P
G03=R3*P
syms s
[Num,Den] = tfdata(G01);
G01_syms = poly2sym(cell2mat(Num),s)/poly2sym(cell2mat(Den),s)
Kv1=limit(s*G01_syms)
[Num,Den] = tfdata(G01);
G01_syms = poly2sym(cell2mat(Num),s)/poly2sym(cell2mat(Den),s)
Kv1=limit(s*G01_syms)
[Num,Den] = tfdata(G02);
G02_syms = poly2sym(cell2mat(Num),s)/poly2sym(cell2mat(Den),s)
Kv2=limit(s*G02_syms)
[Num,Den] = tfdata(G03);
G03_syms = poly2sym(cell2mat(Num),s)/poly2sym(cell2mat(Den),s)
Kv3=limit(s*G03_syms)
%Kv1=10 go ispolnuva samo uslovot
G1=feedback(G01,1)
step(G1)
stepinfo(G1)
% preskokot e 13,14% <20% -> gi ispolnuva proektnite baranja
%% zadaca 3.2
clc 
clear
s=tf('s')
P=1/(s*(s+10)*(s+15))
R1=10/(s+1)
R2=100*((s+1)*(s+5))/((s+10)*(s+50))
R3=2850*(s+1)/(10*s+1)
G01=R1*P
G02=R2*P
G03=R3*P
% A) da ima stac greska pomala od 10%
syms s
[Num,Den] = tfdata(G01);
G01_syms = poly2sym(cell2mat(Num),s)/poly2sym(cell2mat(Den),s)
Kv1=limit(s*G01_syms)
[Num,Den] = tfdata(G01);
G01_syms = poly2sym(cell2mat(Num),s)/poly2sym(cell2mat(Den),s)
Kv1=limit(s*G01_syms)
[Num,Den] = tfdata(G02);
G02_syms = poly2sym(cell2mat(Num),s)/poly2sym(cell2mat(Den),s)
Kv2=limit(s*G02_syms)
[Num,Den] = tfdata(G03);
G03_syms = poly2sym(cell2mat(Num),s)/poly2sym(cell2mat(Den),s)
Kv3=limit(s*G03_syms)
e_inf1=100/Kv1
e_inf2=100/Kv2
e_inf3=100/Kv3 %e3(inf)=5%<10%
%B) ceta=0.707 za dominantniot par polovi
G3=feedback(G03,1)
pole(G3)% dominanten par se -1.2348+-i1.2334
G_tilda_3=1/((s+1.2348-i*1.2334)*(s+1.2348+i*1.2334))
simplify(G_tilda_3)
wn=sqrt(3.046)
syms ceta
simplify(solve(2*ceta*wn==2.46996,ceta))
% ceta=0.7076 sto go ispolnuva uslovot
%% zadaca 3.3
s=tf('s')
P=1/(s*(s/8+1)*(s/20+1))
R1=(100*s+1)/(120*s+1)
R2=((s+1)*(20*s+1))/((100*s+1)*(s/50+1))
R3=100
R4=100*((s+1)*(s/5+1))/((10*s+1)*(s/50+1))
G01=R1*P
G02=R2*P
G03=R3*P
G04=R4*P
% A) Kv>=100
syms s
[Num,Den] = tfdata(G01);
G01_syms = poly2sym(cell2mat(Num),s)/poly2sym(cell2mat(Den),s)
Kv1=limit(s*G01_syms)
[Num,Den] = tfdata(G01);
G01_syms = poly2sym(cell2mat(Num),s)/poly2sym(cell2mat(Den),s)
Kv1=limit(s*G01_syms)
[Num,Den] = tfdata(G02);
G02_syms = poly2sym(cell2mat(Num),s)/poly2sym(cell2mat(Den),s)
Kv2=limit(s*G02_syms)
[Num,Den] = tfdata(G03);
G03_syms = poly2sym(cell2mat(Num),s)/poly2sym(cell2mat(Den),s)
Kv3=limit(s*G03_syms)
[Num,Den] = tfdata(G04);
G04_syms = poly2sym(cell2mat(Num),s)/poly2sym(cell2mat(Den),s)
Kv4=limit(s*G04_syms)
% samo R3 i R4 obezbeduva Kv=100
% B)
G3=feedback(G03,1)
G4=feedback(G04,1)
figure(1)
margin(G03)
figure(2)
margin(G04)

%% zadaca 3.4
s=tf('s')
K=1/0.23
P=(9*K)/((s+1)*(s^2+3*s+9))
nyquist(P)
%% zadaca 3.5
s=tf('s')
K=8.35

P=(K*exp(-0.2*s))/(s+5)
% A)
nyquist(P)
% K=K=1/0.0884
% B)
%% zadaca 3.9
clc
clear
s=tf('s')
P=1/(s+2)
R=10/s
G0=P*R
%margin(G0)

G=feedback(G0,1)

[Peak,freq]=getPeakGain(G)
%% zadaca 3.10
clc
clear
s=tf('s')
G0=100/(s*(s+12))
G=feedback(G0,1)
[Mr,Wr]=getPeakGain(G)
Wo=bandwidth(G)



%% zadaca 3.11
clc
clear
s=tf('s')
syms K
solve(20*log10(K)==-4.08,K)
k= 0.6252
G0=(k*(2*s+sqrt(12)))/s^2
margin(G0)
%% zadaca 3.12
clc
clear
s=tf('s')
G0=(s+3)/(s^2*(s+5))
G=feedback(G0,1)
B=isstable(G)
[Gm,Pm]=margin(G0)
G0c=G0*((s+2)/(s+4))
[Gm1,Pm1]=margin(G0c)
% zadacaa 3.13
G0c1=G0*((s+0.2)/(s+4))
[Gm2,Pm2]=margin(G0c1)
%% zadaca 3.15
s=tf('s')
P=200/((s+2)*(s+4)*(s+5))
Gd1=(s+3)/(s+30)
Gd2=(s+0.3)/(s+30)
G01=Gd1*P
G02=Gd2*P
w0=bandwidth(P)
w01=bandwidth(G01)
w02=bandwidth(G02)

%% zadaca 3.17
clc
clear
s=tf('s')
P=1500/(s*(s+10)*(s+15)*(s+20))
R1=(s+1.2)/(s+0.02)
R2=(s+2)/(s+12)
R3=(2*s+1)/(4*s+1)
R4=(20*s+1)/(120*s+1)
G1=feedback(P*R1,1)
G2=feedback(P*R2,1)
G3=feedback(P*R3,1)
G4=feedback(P*R4,1)
stepinfo(G1)
stepinfo(G2)
stepinfo(G3)
stepinfo(G4)
