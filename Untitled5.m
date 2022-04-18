clc
clear
s=tf('s')
G0=(s+4)/((s+2)*(s+6)*(s+8))

syms ss
G0_syms=(ss+4)/((ss+2)*(ss+6)*(ss+8))
Kp=limit(G0_syms,ss,0)
Kps=100
K=Kps/Kp

K=2400

margin(K*G0)

Phi=-180+45+5

wnew=12.5

b=wnew/10

b_a=10^(21.7/20)

a=b/b_a

I=(a/b)*(s+b)/(s+a)

margin(K*I*G0)
G0_syms=(ss+4)/((ss+2)*(ss+6)*(ss+8))
I_syms=((a/b)*(ss+b)/(ss+a))
Kp=limit(K*I_syms*G0_syms,ss,0)

%% zadaca
clc
clear
clc
clear
s=tf('s')
G0=(s+4)/((s+2)*(s+6)*(s+8))

syms ss
G0_syms=(ss+4)/((ss+2)*(ss+6)*(ss+8))
Kp=limit(G0_syms,ss,0)
Kps=100
K=Kps/Kp

K=2400

margin(K*G0)
Phiwanted=45
Phicurrent=14.1
corr=10

Phimax=Phiwanted-Phicurrent+corr
b_a=(1-sind(Phimax))/(1+sind(Phimax))
20*log(sqrt(1/b_a))
wnew=121

a=wnew*sqrt(1/b_a)
b=wnew*sqrt(b_a)
D=(a/b)*(s+b)/(s+a)

margin(K*D*G0)

%%
clc
clear
s=tf('s')
G0=1/(s*(s+1)*(s+2))
syms ss
G0_syms=1/(ss*(ss+1)*(ss+2))
Kv=limit(ss*G0_syms,ss,0)
Kvs=10
K=Kvs/Kv
K=20
Phi=-180+45+5

