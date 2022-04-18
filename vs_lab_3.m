%% пример 1
a = 5;
b = 0.1;s = tf('s');
Gi1 = (a/b)*((1+s/a)/(1+s/b))
linearSystemAnalyzer(Gi1)

%% пример 2
clc
clear
s=tf('s')
G0=1/(s*(1+s/2)*(1+s/6))
%Kv=5
%w1=1.5rad
%PM e(40,45)

K=5
a=1.34/16
b=a*10^(-9.88/20)
s=tf('s')
G1=K*(1+s/a)/(1+s/b)
margin(G1*G0)

%% пример 3
z=tf('z',0.1)
G0=(3/8)*((z+1)*(z+1/3))/((z-1)*z*(z+1/2))
syms z
[Num,Den] = tfdata(G0);
G0_syms = poly2sym(cell2mat(Num),z)/poly2sym(cell2mat(Den),z)
Kv=limit((z-1)*G0_syms,1)

clc 
clear
z=tf('z',0.1)
G0=(3/8)*((z+1)*(z+1/3))/((z-1)*z*(z+1/2))
G01=7.5*G0

s = tf('s');
G0s = (1/3)*(((1-s)*(1+s/2))/(s*(1+s)*(1+s/3))) 
G0s1 = 7.5*G0s
Gi = (1+s/0.034)/(1+s/0.007) 
G0s2 = Gi*G0s1
ltiview(G0s2)

%% задача 1
s=tf('s')
G0=(5000*(s+4))/((s+1)*(s^2+16*s+100)*(s+20))
%syms s
%[Num,Den] = tfdata(G0);
%G0_syms = poly2sym(cell2mat(Num),s)/poly2sym(cell2mat(Den),s)
%Kp=limit(G0_syms)
K=4
G01=K*G0
%ltiview(G01)
% Od najkvitovata kriva se gleda deka sistemot e nestabilen
% i mora da se stabilizira so integralen komp
%od toa sto se bara w1=15, kje se zeme nulata 10 pati pomala od w1
a=15/10
b=a*10^(-9.71/20)
% za frek od 15 zasiluvanjeto e 9.71,pa znaci deka treba 
%da se vnese zasiluvanje od -9.71
Gi=(1+s/a)/(1+s/b)
G02=Gi*G01
ltiview(G02)
% znaci so voveduvanje na vakvite promeni,se kompenzira sistemot da gi
% ispolnuva proektnite baranja

%% ЗАДАЧА 2
s=tf('s')
G0=(s+4)/((s+2)*(s+6)*(s+8))
%syms s
%[Num,Den] = tfdata(G0);
%G0_syms = poly2sym(cell2mat(Num),s)/poly2sym(cell2mat(Den),s)
%Kp=limit(G0_syms)
%K=100/Kp 
K = 2400
G01=K*G0
%ltiview(G01)
% od najkvist se gleda deka vakviot sistem e stabilen, ostanuva da se
% proveri dali go ispolnuva drugoto proektno baranje
%od bode jasno e deka ne gi ispolnuva, treba da se kompenzira so soodveten
%integralen komp
a=0.2
b=0.02
Gi=(1+s/a)/(1+s/b)
G02=Gi*G01
ltiview(G02)
% so proba i greska stigam do a=0.2 i b=0.02, pri sto se ovozmovuva
% PM=44.3~45,ostanuva da proveram dali e stabilen
% sistemot e stabilen

%% задача 3
clc
clear
z=tf('z',1)
G0z=(27/64)*((z+1)/(z+0.5))^3
s=tf('s')
G0=d2c(G0z,'tustin')
%go prefrlam sistemot vo kontinualen domen so bilinearna transformacija
%i kako takov go proektiram
%syms s
%[Num,Den] = tfdata(G0);
%G0_syms = poly2sym(cell2mat(Num),s)/poly2sym(cell2mat(Den),s)
%Kp=limit(G0_syms)
Kp=1
K=4
G01=4*G0
%ltiview(G01)
a=3.8/10
b=a*10^(-7.59/20)
Gi=(1+s/a)/(1+s/b)
G02=Gi*G01
ltiview(G02)
% so proba i greska stigam do PM=78.7 I GM=13.3 so sto se ispolneti
% proektnite baranja
%% задача 4
clc
clear
s=tf('s')
G0=120*((s+1)/(s*(s^2+s+4)*(s+8)*(s+16)))
%syms s
%[Num,Den] = tfdata(G0);
%G0_syms = poly2sym(cell2mat(Num),s)/poly2sym(cell2mat(Den),s)
%Kv=limit(s*G0_syms) 
%K=30/Kv
K=128

G01=K*G0
%ltiview(G01)
% od najkvist se gleda deka sistemot e nestabilen pa mora da se kompenzira
% nezavisno od toa dali gi ispolnuva baranjata
a=3/20


b=a*10^(-29/20)
Gi=(1+s/a)/(1+s/b)
G02=Gi*G01
%ltiview(G02)

c=3
d=33
Gd=(1+s/c)/(1+s/d)
G03=Gd*G02
ltiview(G03)
