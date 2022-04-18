%% zadaca 2
clc
clear
s=tf('s')
G0=1/(10000*(s^2-1.1772))
rlocus(G0)
sgrid(0.7,0.5)

theta_poz=180-atand(0.357/(1.1772+0.35))
theta_neg=atand(0.357/(1.1772-0.35))
psi_a=theta_poz+theta_neg-180

a=0.357/tand(psi_a)+0.35
Kd=5.78*10^3
D=Kd*(s+a)

rlocus(D*G0)
sgrid(0.7,0.5)