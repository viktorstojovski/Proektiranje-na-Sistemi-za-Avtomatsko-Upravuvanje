%% zadaca 3
clc
clear
s=tf('s')
Mp=30
ts=5
zeta=-log(Mp/100)/sqrt(pi^2+(log(Mp/100))^2)
wn=4/(zeta*ts)

P=1/(s*(s^2+2))
G=feedback(P,1)
syms s
as=expand((s+zeta*wn+i*wn*sqrt(1-zeta^2))*(s+zeta*wn-i*wn*sqrt(1-zeta^2))*(s+5*zeta*wn))
as=s^3+5.6*s^2+11.4*s+20

Aw=[0 1 0; 0 0 1;-1 -2 0]
Bw=[0;0;1]
Cw=[1 0 0]

syms k1 k2 k3
K=[k1 k2 k3]
K=[19 9.4 5.6]
det(s*eye(3)-Aw+Bw*K)

Av=transpose(Aw)
Bv=transpose(Cw)
Cv=transpose(Bw)
als=expand((s+8+i*20.875)*(s+8-i*20.875)*(s+40))
als=s^3+56*s^2+1139.8*s+19991
syms lv1 lv2 lv3
Lv=[lv1;lv2;lv3]
det(s*eye(3)-Av+Lv*Cv)
Lv=[19990;1137.8;56]

Qov=obsv(Av,Cv)
Qow=obsv(Aw,Cw)

P=Qov*inv(Qow)

Lw=P*Lv
eig(Aw-Bw*K)
eig(Aw-Lw*Cw)