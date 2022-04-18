%% zadaca 1
clc
clear
s=tf('s')
G0=40/(s*(s+2))
margin(G0)

wnew=10
Phiwanted=30
Phicurrent=18
corr=1
Phimax=Phiwanted-Phicurrent+corr

b_a=(1-sind(Phimax))/(1+sind(Phimax))
20*log(sqrt(1/b_a))
a=wnew*sqrt(1/b_a)
b=wnew*sqrt(b_a)

D=(a/b)*(s+b)/(s+a)

margin(D*G0)
