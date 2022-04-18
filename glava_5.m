
%% zadaca 5.12
clc
clear
s=tf('s')
K=10^(3.64/20)
P=K*exp(-2*s)/(s+1)
hold on
margin(P)
%step(feedback(P,1))
%% zadaca 5.13
clc
clear
s=tf('s')
P=1/((s+1)*(s+3)*(s+5))
K1=10^(-8.3/20)
K=500*K1


G0=K*P
margin(G0)
%% zadaca 5.16
clc
clear
s=tf('s')
P=1/(s*(s+10))
K=10^(43.1/20)
G0=K*P
margin(G0)
%% zadaca 5.17
clc
clear
s=tf('s')
P=1/(s*(s+5)^2)
K=10^(35.5/20)
G0=K*P
margin(G0)
%% zadaca 5.21
clc
clear
s=tf('s')
P=10/(s*(s+1))
R=25/(s+10)

K=10^(-28/20)
G0=K*P*R
margin(G0)
figure
step(feedback(G0,1))
%% zadaca 5.22
clc
clear
s=tf('s')
P=1/(s*(s+10))
K=200


%D=(b/a)*(s+a)/(s+b)
Phicurrent=38.7
Phiwanted=45
corr=4
Phimax=Phiwanted-Phicurrent+corr
a_b=(1-sind(Phimax))/(1+sind(Phimax))%a/b

20*log10(sqrt(1/a_b))
wnew=14.1
a=wnew*sqrt(a_b)
b=wnew*sqrt(1/a_b)
D=(b/a)*(s+a)/(s+b)
G0=K*P*D
margin(G0)

%% zadaca 5.23
clc
clear
s=tf('s')
G0=72/((s+1)*(s+3)^2)
% D=(sqrt(b/a))*(s+a)/(s+b)

wnew=3.38
Phiwanted=45
Phicurrent=9.66
corr=0

Phimax=Phiwanted-Phicurrent+corr

a_b=(1-sind(Phimax))/(1+sind(Phimax))
a=wnew*sqrt(a_b)
b=wnew*sqrt(1/a_b)
 D=(sqrt(b/a))*(s+a)/(s+b)
Gc=G0*D
figure(1)
margin(G0)
figure(2)
margin(D)
figure(3)
margin(Gc)
figure(4)
step(feedback(Gc,1))

%% zadaca 5.24
clc
clear
s=tf('s')
G0=100/(s*(s+50)*(s+100))

K=2000
G0_1=K*G0

margin(G0_1)

Phiwanted=48
Phicurrent=39.6
corr=10

Phimax=Phiwanted-Phicurrent+corr
a_b=(1-sind(Phimax))/(1+sind(Phimax))%a/b

20*log10(sqrt(1/a_b)) %barame pri koja frekv , zasiluvanjeto e  -2.8386

wnew=40%taa frekv e 40
a=wnew*sqrt(a_b)
b=wnew*sqrt(1/a_b)
D=(b/a)*(s+a)/(s+b)
G0_2=G0_1*D
figure(1)
%margin(G0_2)
figure(2)
%step(feedback(G0_2,1))
syms s
[Num,Den] = tfdata(G0_2);
G01_syms = poly2sym(cell2mat(Num),s)/poly2sym(cell2mat(Den),s)
Kv1=limit(s*G01_syms)% =40
%% zadaca 5.25
clc
clear
s=tf('s')


G0=32/(s*(s+2)*(s+8))

% nagoduvanje na K

syms ss
G0_syms=32/(ss*(ss+2)*(ss+8))
Kv=limit(ss*G0_syms,ss,0)
e_inf=(10*pi)/Kv
K=50
figure(1)
margin(G0)

Phiwanted=45
Phicurrent=-43
corr=0

Phimax=Phiwanted-Phicurrent+corr

a_b=(1-sind(Phimax))/(1+sind(Phimax))

20*log10(sqrt(1/a_b))% =35.16

wnew=w0





%% zadaca 5.26
clc
clear
s=tf('s')
G0=5/(s*(s+1)*(s+2)*(s+3))
figure(1)
margin(G0)

Phiwanted=50
Phicurrent=26.8
corr=10


Phimax=Phiwanted-Phicurrent+corr
az
a_b=(1-sind(Phimax))/(1+sind(Phimax))

20*log10(sqrt(1/a_b))

wnew=0.95
a=wnew*sqrt(a_b)
b=wnew*sqrt(1/a_b)
D=(b/a)*(s+a)/(s+b)

 
figure(2)
margin(D*G0)



%% zadaca 5.27
clc
clear
s=tf('s')
G0=100/(s*(s+50)*(s+100))
Mp=10
epsi=(-log(Mp/100))/(sqrt(pi^2+((-log(Mp/100)))^2))
PM=atand(2*epsi/sqrt(-2*epsi^2+sqrt(1+4*epsi^4))) %58.59

PMmax_K=-180+PM %=-121.406
figure(1)
margin(G0)
%so cel da se zadrzi ist Mp,pri koj PM treba da e 58.59,potrebno e
%zasiluvanjeto da bide 
K1=10^(3)

syms ss
G0_syms=K1*100/(ss*(ss+50)*(ss+100))
Kv=limit(ss*G0_syms,ss,0)
e_inf=1/Kv
e_s_inf=e_inf/10
K=10

G0_K_komp=K*K1*G0

figure(2)
margin(G0_K_komp)
wnew=18.5% magnitudata na zasiluvanjeto e 20

b=wnew/10
b_a=10^(20/20)
a=b/b_a

I=(a/b)*(s+b)/(s+a)

figure(3)
margin(I*G0_K_komp)

%% zadaca 5.28
clc
clear
s=tf('s')
G_nc=1/(s*(s+20)^2)
epsi=0.707
PM=atand(2*epsi/sqrt(-2*epsi^2+sqrt(1+4*epsi^4)))% 65.5246
syms s
[Num,Den] = tfdata(G_nc);
G01_syms = poly2sym(cell2mat(Num),s)/poly2sym(cell2mat(Den),s)
Kv1=limit(s*G01_syms)
s=tf('s')
K=40*400
G_Kc_Inc=K*G_nc
margin(G_Kc_Inc)
Phimax=-180+65.5+5% -109.5 pri 3.58, a mag e 20.7dB
wnew=3.58
a=wnew/10
a_b=10^(20.7/20)
b=a/a_b
I=(b/a)*((s+a)/(s+b))
G_c=I*G_Kc_Inc
%margin(G_c)

%% zadaca 5.29
clc
clear
s=tf('s')
G0=1/(s*(s+5)^2)

Kv1=1/25

K=5*(1/Kv1)+1
%margin(K*G0)

wnew=1.01
a=wnew/10
a_b=10^(13.6/20)
b=a/a_b
I=(b/a)*((s+a)/(s+b))
margin(I*K*G0)

%% primer 5.11
clc
clear
s=tf('s')
K=1

K=20/4
P=(4*K)/(s*(s+1))

syms ss
P_syms=4/(ss*(ss+1))
Kv=limit(ss*P_syms,ss,0)


K=20/4


margin(P)
Phiwanted=35
corr=5
Phimax=-180+Phiwanted+corr % iznesuva -140 pri frek od 1.2

wnew=1.1 % zasiluvanjeto e 20.6

a_b=10^(20.6/20)

a=wnew/10
b=a/a_b

I=(b/a)*((s+a)/(s+b))
margin(I*P)
%% primer 5.15
clc
clear
s=tf('s')
G0=8/(s*(s+3)*(s+5))

%proektni baranja
% 1) Mp<=10%
% 2)Tm<=1s
% 3) Kv>=20

% nagoduvanje na K

syms ss
G0_syms=8/(ss*(ss+3)*(ss+5))
Kv=limit(ss*G0_syms,ss,0)

Kvs=20

K=Kvs/Kv
K=37.5
% iscrtuvanje na K*G0
margin(K*G0)


% presmetuvanje na vrski na preoden rezim
Mp=10
zeta=(-log(Mp/100))/(sqrt(pi^2+((-log(Mp/100)))^2))
Tm=1
w0=(pi/(Tm*sqrt(1-zeta^2)))/sqrt((1-2*zeta^2+sqrt(4*zeta^4-4*zeta^2+2)))

PM=atand(2*zeta/sqrt(-2*zeta^2+sqrt(1+4*zeta^4)))
corr=5
Phimax=PM-180+172+corr

wnew=w0
b2=wnew/10 %nula
b_a=(1+sind(Phimax))/(1-sind(Phimax))

a2=b2/b_a

b1=wnew*sqrt(1/b_a)
a1=wnew*sqrt(b_a)

I=(a2/b2)*(s+b2)/(s+a2)
D=(a1/b1)*(s+b1)/(s+a1)
figure(1)
margin(K*I*D*G0)
figure(2)
step(feedback(K*I*D*G0,1))

syms s
[Num,Den] = tfdata(K*I*D*G0);
G01_syms = poly2sym(cell2mat(Num),s)/poly2sym(cell2mat(Den),s)
Kv1=limit(s*G01_syms)

%% zadaca 5.30

clc
clear
s=tf('s')
G0=1/(s*(s+8)*(s+30))

%proektni baranja
% 1) Mp=10%
% 2)Tp=0.6
% 3) Kv=10


%nagoduvanje na Kv

syms ss
G0_syms=1/(ss*(ss+8)*(ss+30))
Kv=limit(ss*G0_syms,ss,0)
Kvs=10
K=Kvs/Kv
K=2400

%iscrtuvanje na bode za K*G0
margin(K*G0)

% presmetuvanje vrski

Mp=10
zeta=(-log(Mp/100))/(sqrt(pi^2+((-log(Mp/100)))^2))
Tm=0.6
w0=(pi/(Tm*sqrt(1-zeta^2)))/sqrt((1-2*zeta^2+sqrt(4*zeta^4-4*zeta^2+2)))
PM=atand(2*zeta/sqrt(-2*zeta^2+sqrt(1+4*zeta^4)))
wnew=w0


Phimax=PM-180+129+12

b_a=(1+sind(Phimax))/(1-sind(Phimax))
b2=wnew/10
a2=b2/b_a

I=(a2/b2)*(s+b2)/(s+a2)

b1=wnew*sqrt(1/b_a)
a1=wnew*sqrt(b_a)

D=(a1/b1)*(s+b1)/(s+a1)


figure(1)
margin(K*I*D*G0)
figure(2)
step(feedback(K*I*D*G0,1))

%% zadaca (1 zadaca od lab4)
clc
clear
s=tf('s')


G0=25/((s+0.14)*(s+4)*(s+15))

%proektni baranja
%1)Kp=100
%2)w1=5rad/s
%3)PM=30

% nagoduvanje na stacionaren rezim

syms ss
G0_syms=25/((ss+0.14)*(ss+4)*(ss+15))
Kp=limit(G0_syms,ss,0)
Kps=100
K=Kps/Kp
K=33.6

%iscrtuvanje na Bode za K*G0
figure(1)
margin(K*G0)

wnew=5

Phimax=30-180+158+10

b_a=(1+sind(Phimax))/(1-sind(Phimax))
b2=wnew/10
a2=b2/b_a

I=(a2/b2)*(s+b2)/(s+a2)

b1=wnew*sqrt(1/b_a)
a1=wnew*sqrt(b_a)

D=(a1/b1)*(s+b1)/(s+a1)
figure(2)
margin(K*I*D*G0)



%% zadaca (4 od lab 4)
clc
clear
s=tf('s')
G0=120*(s+1)/(s*(s^2+s+4)*(s+8)*(s+16))

%proektni baranja:
%1)Kv=30
%2)Pm>=25
%3)w1=3rad/s
%4)wpi>=8

% nagoduvanje na stacionaren rezim

syms ss
G0_syms=(120*(ss+1))/(ss*(ss^2+ss+4)*(ss+8)*(ss+16))
Kv=limit(ss*G0_syms,ss,0)
Kvs=30
K=Kvs/Kv
K=128

%iscrtuvanje na bode za K*G0
figure(1)
margin(K*G0)

% Dizajniranje na I kompenzator
wnew=3
%presmetuvanje na Phimax
Phimax=25-180+228+10

%presmetuvanje na odnos nula pol

b_a=(1+sind(Phimax))/(1-sind(Phimax))
b2=wnew/10
a2=b2/b_a

I=(a2/b2)*(s+b2)/(s+a2)

%Dizajniranje na D kompenzator

b1=wnew*sqrt(1/b_a)
a1=wnew*sqrt(b_a)

D=(a1/b1)*(s+b1)/(s+a1)

%iscrtuvanje na bode za kompenziraniot sistem
figure(2)
margin(K*I*D*G0)

%proverka da ne se slucilo neso vo stac rezim


syms s
[Num,Den] = tfdata(K*I*D*G0);
G01_syms = poly2sym(cell2mat(Num),s)/poly2sym(cell2mat(Den),s)
Kv1=limit(s*G01_syms)




%% zadaca (5 od lab 4)
clc
clear

