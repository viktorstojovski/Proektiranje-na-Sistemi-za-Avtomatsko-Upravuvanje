%% Zadaca 5.12
clc
clear
s=tf('s')
K=500*10^(-8.31/20)
P=1/((s+1)*(s+3)*(s+5))
margin(K*P)