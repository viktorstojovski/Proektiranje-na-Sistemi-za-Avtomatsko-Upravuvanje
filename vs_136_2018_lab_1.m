%% zadaca 1
%proektni baranja:
%vreme na smiruvanje <2s
%preskok<5%
%stacionarna greska <1% pri ref vlez od 1 rad
J = 3.2284E-6;
b = 3.5077E-6;
K = 0.0274;
R = 4;
L = 2.75E-6;
s = tf('s');
P_motor = K/(s*((J*s+b)*(L*s+R)+K^2));
Kp=13

Ki=0
Kd=0
C = pid(Kp,Ki,Kd)
H=feedback(C*P_motor,1)
stepinfo(H)
step(H)

syms s
[Num,Den] = tfdata(C*P_motor);
G0_syms = poly2sym(cell2mat(Num),s)/poly2sym(cell2mat(Den),s)
Kp=limit(G0_syms)
e=1/Kp
%% zadaca 2
clc
clear
M=1
b=10
k=20
s = tf('s')
P=1/(M*s^2+b*s+k)
Kp=0.0040079

Ki=0.001235
Kd=0.00030875
C = pid(Kp,Ki,Kd)
G=feedback(C*P,1)
stepinfo(G)
syms s
[Num,Den] = tfdata(C*P);
G0_syms = poly2sym(cell2mat(Num),s)/poly2sym(cell2mat(Den),s)
Kp=limit(G0_syms)
e=1/Kp

%% zadaca 3
clc
clear
M=1000
b=50

s = tf('s');
P=1/(M*s+b)

Kp=6000


Ki=0
Kd=10

C = pid(Kp,Ki,Kd)
G=feedback(C*P,1)
stepinfo(G)
syms s
[Num,Den] = tfdata(C*P);
G0_syms = poly2sym(cell2mat(Num),s)/poly2sym(cell2mat(Den),s)
Kp=limit(G0_syms)
e=2/Kp

%% zadaca 4
clc
clear
s = tf('s');
G=0.5/(s^2+s+1)
Kp=100


Ki=10
Kd=0

C = pid(Kp,Ki,Kd)
Gn=G/(1-G*C)
En=Gn/s
syms s
syms s
[Num,Den] = tfdata(En);
En_syms = poly2sym(cell2mat(Num),s)/poly2sym(cell2mat(Den),s)
en=limit(s*En_syms,s,0)
%moze da se zabelezi deka ne vazno od parametrite na pi upravuvacot
%,greskata kje e nuleva ako se upotrebi pi upravuvac

%b
E1n=Gn-1
E=En+E1n
% ne razbiram

