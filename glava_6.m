%% zadaca 6.28
clc
clear
%iscrtuvanje na GMK
s=tf('s')
G0=1/((s+1)*(s+15)^2)
rlocus(G0)
%presmetka na Kp
syms s
G0=1/((s+1)*(s+15)^2)
Kp=limit(G0,s,0)
% Kp=K/225>=20 -> K=Kp*20>=4500 granica (1)

%presmetuvanje na d=Kgranicno/Kkonstructivno
Kgranicno=7.67*10^3 %se opredeluva od GMK
%d=Kgranicno/Kkonstructivno>=1.5 -> K=Kkonstruktivno=Kgranicno/1.5
K=Kgranicno/1.5 %=5.1133e+03
% 4500<=K<=5113.3

%% zadaca 6.30
clc
clear
s=tf('s')
G0=20/(s*(s+10)^2)
rlocus(G0)
a=0.0001
b=100*a
I=(s+b)/(s+a)
hold on
rlocus(I*G0)

syms s
G0=20/(s*(s+10)^2)
Kv=limit(s*G0,s,0)
%  Kv =20/100
%lamda=100 e odnosot Kvs/Kv=100, odnosno es(inf)/e(inf)=1/100, a toa gi
%definira i odnosot na polot i nulata na integralniot b/a=100
%tie treba da se mnogu blisku do imaginarnata oska najdesno od site nuli i
%polovi isklucuvajki go polot vo 0

%% zadaca 6.32
%vo ovoj primer potrebno e da se nagodi i preodniot i stacionarniot rezim
%redosledot po koj toa se vrsi e 
% 1)preoden rezim (naogjanje na K pri koj dom par polovi kje ima zeta=0.592
% 2)stacionaren rezim (naogjannje na I kompenzator
clc
clear
s=tf('s')
G0=1/(s*(s+5))

%1) nagoduvanje na preoden rezim
rlocus(G0)

% K za zeta=0.592 e K=17.8, a dom.par polovi se s1/2=-2.5±j3.39
K=17.8
% 2) nagoduvanje na stacionaren rezim
syms s
G0_syms=17.8/(s*(s+5))
Kv=limit(s*G0_syms,s,0)
e_inf_nc=1/Kv
% opredeluvanje na odnosot b/a=lambda

%lambda=Kvs/Kv
%za e_inf_c=0.02 Kvs=1/0.02
Kvs=1/0.02
b_a=Kvs/Kv
% se dobiva 14.0449~14
b_a=14
b=0.1
a=b/b_a
s=tf('s')

I=(s+b)/(s+a)
hold on
rlocus(K*I*G0)
sgrid(0.592,[])


%% zadaca 6.33
clc
clear

s=tf('s')
K=3
G0=K/(s*(s+2)^2)

rlocus(G0)
syms s
G0_syms=3/(s*(s+2)^2)
Kv=limit(s*G0_syms,s,0)

% Kv=3/4 e_inf=4/3=1.3333
%so cel da se namali greskata 50 pati potrebno e da se proektira
%realen kompenzator so odnos b_a=lambda=50 , no najprvo daj da go iscrtame
%GMK, nulata i polot treba da se megju(-0.51,0)
%neka b=0.1 , a=b/b_a
s=tf('s')
b=0.1
a=b/50
I=(s+b)/(s+a)

hold on
rlocus(I*G0)
syms s
I_G0_syms =3*(s+0.1)/((s*(s+2)^2)*(s+a))
Kv=limit(s*I_G0_syms,s,0)


%% zadaca 6.34
%dokolku e potrebna nuleva stac. greska treba da se iskoristi idealen I
%kompenzator od oblik I=(s+b)/s
clc
clear
s=tf('s')
G0=1/((s+2)*(s+4)*(s+8))
%1) nagoduvanje na preoden rezim
%rlocus(G0)
sgrid(0.18,[]) 
% za zeta=0.18, K=316 i s1/2=-1.01±j5.54
K=316
%2) nagoduvanje na stac
%nulata se bira blisku do Im{s}, na primer b=-0.005

I=(s+0.005)/s

hold on
rlocus(K*I*G0)

%% zadaca 6.35
clc
clear

clc
clear
s=tf('s')
G0=1/((s+2)*(s+4)*(s+8))



%syms s
%G0_syms=314/((s+2)*(s+4)*(s+8))
%Kp=limit(G0_syms,s,0) %
Kp=  4.9062
e_inf_nc=1/(1+Kp) %e_inf_nc=0.1693
e_inf_c=e_inf_nc/10
Kp_c=1/e_inf_c-1
b_a=Kp_c/Kp
b=0.118
a=b/b_a

s=tf('s')
I=(s+b)/(s+a)

rlocus(314*I*G0)
step(feedback(314*I*G0,1))
stepinfo(feedback(314*I*G0,1))



%% zadaca 6.36
clc
clear
s=tf('s')
G0=1/(s*(s+2))

zeta=-log(53.3/100)/(sqrt(pi^2+(log(53.3/100))^2))% =  0.1964
rlocus(G0)
sgrid(0.2044,[])% pri zeta=0.2044,Zasiluvanjeto K=24
K=24
%a) vreme na smiruvanje?
stepinfo(feedback(K*G0,1))
Ts=3.4988 %sec

%b)stac. greska pri edinecno rastecki vlez
syms ss
G0_syms=K/(ss*(ss+2))
Kv=limit(ss*G0_syms,ss,0)
e_inf=1/Kv %=1/12

%v)I-komp koj kje ostvari e_inf_c=0.01
e_inf_c=0.01
e_inf_nc_c=e_inf/e_inf_c
Kvs=e_inf_nc_c*Kv

b_a=Kvs/Kv
b=0.1
a=b/b_a

I=(s+0.1)/(s+0.0120)

G0_c=K*I*G0

step(feedback(G0_c,1))

G_c=feedback(G0_c,1)
zpk(G_c)

G_c_2_red=24/(s^2 + 1.911*s + 23.83)
hold on
step(G_c_2_red)
%% zadaca 6.37
clc
clear
s=tf('s')
G0=1/(s*(s+3)*(s+5))
stepinfo(feedback(G0,1))
% opredeluvanje na K, i dominantnite polovi
zeta=-log(16/100)/(sqrt(pi^2+(log(16/100))^2))
rlocus(G0)
sgrid(zeta,[])

% od GMK moze da se otcita deka K=21.3, a dominantniot par na polovi
% s1/2=-0.941±j1.61
Ts=4/0.941 %Ts~4sec

Tsc=Ts/4
Tsc=1


% za da se postigne toa potrebno e dom. par polovi na zatvoreniot sistem da
% se ednakvi na sc1/2=4*s1/2=-3.764±j 6.44

theta_0=180-atand(6.44/3.764)
theta_3=180-atand(6.44/0.764)
theta_5= atand(6.44/(5-3.764))


psi_a=-(180-theta_0-theta_3-theta_5)

% a: tan(psi_a)=6.44/(3.764-a); a=3.764-6.44/tan(.)
a=3.767-6.44/abs(tand(psi_a))
Kd=40
D=Kd*(s+a)
rlocus(D*G0)
zpk(feedback(D*G0,1))

GC=40/(s^2+7.536*s+51.5)
step(GC)

%% zadaca 6.38
clc
clear
s=tf('s')
G0=1/(s*(s+0.8))

%a) Ts-?

rlocus(G0)
ylim([-2000 2000])

%za K=64 s1/2=-0.4±j8
Ts=4/0.4

%b) ev_inf-?
syms ss
G0_syms=64/(ss*(ss+0.8))
Kv=limit(ss*G0_syms,ss,0)

%v) da se proektira D upravuvac za Ts da se namali 5 patii da se zeme
%nulata na kompenzatorot vo s=-1

a=1

%za da se namali Ts 5 pati, potrebno e dominantniot par na polovi da se
%zgolemi 5 pati odnosno s_c1/2=-2±j40

theta_0=180-atand(20/2)
theta_08=180-atand(20/1.2)
psi_1=180-atand(20)

theta_b=180+psi_1-theta_0-theta_08
b=20/tand(theta_b)+2
Kd=404

D=Kd*((s+a)/(s+b))

rlocus(D*G0)
ylim([0 30])

% iskaca deka Kd=404
%% zadaca 6.40
clc
clear
s=tf('s')
G0=1/(s*(s+5))

zeta=-log(20/100)/(sqrt(pi^2+(log(20/100))^2)) % =0.5965 s1/2=-2.5±j3.336

rlocus(G0)
sgrid(0.5965,[])
K=17.6
%a) Ts-?
stepinfo(feedback(K*G0,1))
Ts=1.4148  %presmetano od dom par polovi =1.6
Ts=1.6
%b) ev(inf)=?
syms ss
G0_syms=17.6/1/(ss*(ss+5))
Kv=limit(ss*G0_syms,ss,0)

%v) D->Tsc=Ts/2; I->ec(inf)=ev(inf)/10

% proektiranje na D upravuvac
Tsc=Ts/2

% s1/2_c=2*s1/2=-5±j6.72

% princip na modul neka ja postavime nulata vo a=-6
theta_0=180-atand(6.72/5)
theta_5=90
psi_a=atand(6.72)
theta_b=180-(theta_0+theta_5-psi_a) %theta_b=44.885

%tan(theta_b)=6.72/(b-5)->
b=6.75/tand(theta_b)+5
% b~-12
Kd=70

Gd=Kd*((s+6)/(s+12))
 
rlocus(Gd*G0)
%) proektiranje na I upravuvac
% b_a=10
b_a=10
b=0.01
a=b/b_a

I=(s+b)/(s+a)

%proverka
Gc0=Gd*I*G0
G=feedback(Gc0,1)
zpk(G)
%se gleda deka Ts e namaleno 2 pati
stepinfo(G)

G_D_c_syms=(Kd*((ss+6)/(ss+12)))*G0_syms
Kvci=limit(ss*G_D_c_syms,ss,0)
ec_inf=1/Kvci

G_c_syms=((ss+b)/(ss+a))*(Kd*((ss+6)/(ss+12)))*G0_syms
Kvc=limit(ss*G_c_syms,ss,0)
ec_inf=1/Kvc

Kvc/Kvci==10 %true





