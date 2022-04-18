%% zadaca 1
clc
clear
s = tf('s');

G=(20)/(s*(1+s/25))
margin(G)
syms s
[Num,Den] = tfdata(G);
G0_syms = poly2sym(cell2mat(Num),s)/poly2sym(cell2mat(Den),s)
Kv=limit(s*G0_syms)
% sistemot nema potreba da se kompenzira, bidejki gi ispolnuva baranjata
% samiot po sebe

%% zadaca 2
clc
clear
s = tf('s');
G=(4)/(s*(s+2))
Gd=100*((s+1.57)/(s+10*1.57))
syms s
margin(Gd*G)
[Num,Den] = tfdata(Gd*G);
G0_syms = poly2sym(cell2mat(Num),s)/poly2sym(cell2mat(Den),s)
Kv=limit(s*G0_syms)
%% zadaca 3

