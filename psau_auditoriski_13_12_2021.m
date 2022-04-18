%% Primer 1
clc
clear

s = tf('s');

G0 = (s+0.1)/(s*(s+0.2))

plot(-0.1, 0, 'o', 'color', 'blue')
hold on;
plot(-0.2, 0, 'x', 'color', 'red')
hold on;
plot(0, 0, 'x', 'color', 'red')
grid on;
xlim([-0.5, 0.2])
ylim([-0.5, 0.5])

% broj na asimptoti = n - m = 2 - 1 = 1
% centar => sigma_c = (0.2 - 0.1)/1 = 0.1
% B1 => (2*0 + 1)*Pi / 1 => Pi

rlocus(G0)

%% Primer 2
clc
clear

syms s
G0(s) = (s+1)/((s+2+3j)*(s+2-3j))

plot(-1, 0, 'o', 'color', 'blue')
hold on;
plot(-2, 3, 'x', 'color', 'red')
hold on;
plot(-2, -3, 'x', 'color', 'red')
grid on;
xlim([-6, 1])
ylim([-4, 4])

% broj na asimptoti = n - m = 2 - 1 = 1
% centar => sigma_c = (2 + 3j + 2 - 3j - 1)/1 = 3
% B1 => (2*0 + 1)*Pi / 1 = Pi

% tochka na razdv. :

syms sigma_b

sol = double(solve(1/(sigma_b + 2+3j) + 1/(sigma_b + 2 - 3j) - 1/(sigma_b + 1) == 0))

plot(sol(1), 0, '*', 'color', 'cyan')

% agol na oddalechuvanje

G0_p1(s) = (s+1)/(s+2+3j)
G0_p2(s) = (s+1)/(s+2-3j)

Theta_d1 = vpa(180 + rad2deg(angle(G0_p1(-2+3j))), 5)
Theta_d2 = vpa(180 + rad2deg(angle(G0_p2(-2-3j))), 5)

s = tf('s');

G0 = (s+1)/((s+2+3j)*(s+2-3j))

rlocus(G0)
hold on;

plot(-2:0.01:10, 3*ones(size(-2:0.01:10, 2), 1), 'k--')
hold on;
plot(-2:0.01:10, -3*ones(size(-2:0.01:10, 2), 1), 'k--')

hold on;

th = 0:1/360:deg2rad(Theta_d1);  r = 0.3;
xx = r*cos(th); yy = r*sin(th);
plot(-2 + xx, 3 + yy, 'r')

th = 0:1/360:deg2rad(Theta_d2);  r = 0.3;
xx = r*cos(th); yy = r*sin(th);
plot(-2 + xx, -3 + yy, 'r')

%% 6.29
clc
clear

s = tf('s');

K = 2;
G0 = K * (s+2)/(s^2 + 4*s + 5)

pole(feedback(G0, 1))

%% 6.28
clc
clear

% a)
syms s K

G0 = K/((s+1)*(s+15)^2)

Kp = limit(G0, s, 0)

sol = solve(Kp >= 20, 'ReturnConditions', true)

sol.conditions

% za Kp >= 20 potrebno e K >= 4500 (uslov 1)
clc

s = tf('s');
G0 = 1/((s+1)*(s+15)^2)

% rezervrata na zasiluvanje d se definira kako: d = K_gr/K_kontruktivno

% K_gr e vrednosta za K za koja sistemot se naogja na granica na stabilnost
% K_konstr. e potrebnata vrednost na K za da se dobie baranata r.z. d

% d = K_gr/K_kontruktivno => K_konstr. = K_gr/d

% za da se najde K_gr se koristi rlocus
rlocus(G0)

% od rlocus zakluchuvame deka sistemot se naogja na gr. na stab. pri
% K_gr ~= 7680

d = 1.5
K_gr = 7680

syms K_konstr

sol = solve(K_gr/K_konstr >= 1.5, 'ReturnConditions', true)

sol.conditions

% za da d >= 1.5 potrebno e K <= 5120 (uslov 2)

% sleduva deka 4500 <= K <= 5120

%% 6.25
clc
clear

s = tf('s');

G0 = 1/(s*(s^2 + 2*s + 5))

rlocus(G0)

% od rlocus zabelezhuvame deka sistemot se naogja na gr. na stab.
% pri a = 10

%% 6.23
clc
clear

s = tf('s');

G0 = (s^2 + 10*s + 74)/(s^2 * (s+4))

rlocus(G0)
sgrid(0.512, 50)

% ceta = 0.512 pri K = 49.2 i se dobivaat dominantni polovi
% s_1/2 = -4.66 +- j7.82

K = 49.2;
pole(feedback(G0 * K, 1))

%% 6.24
clc
clear

syms K s G0

G = K/((s+25)^2 + K);

G0 = solve(G == G0/(1+G0), G0);

% G0 se dobiva deka e G0 = K/(s+25)^2

s = tf('s');

G0 = 1/(s+25)^2

rlocus(G0)
sgrid(0.707, 50)
ylim([-30 30])

% za ceta = 0.707 K treba da bide K = 625 i se dobivaat d.p.
% s_1/2 = -25 +- j25

K = 625;
pole(feedback(G0*K, 1))

%% 6.30 - neidealen I kompenzator
clc
clear

syms s

G0 = 20/(s*(s+10)^2)

Kv = double(limit(s*G0, s, 0))

e_stac = 1/Kv

e_stac_new = e_stac / 100

Kv_new = 1/e_stac_new

gain_needed = double(Kv_new/Kv) % od ova go naogjame potrebniot soodnost z_p na komp.

% Gi(s) = (s+z)/(s+p), z > p

z_p = gain_needed;

z = 0.1; % proizvolno, vo blizina na imaginarnata oska
p = z * 1/z_p % se presmetuva spored potrebniot soodnos i proizvolnata nula

s = tf('s');

G0 = 20/(s*(s+10)^2)
Gi = (s+z)/(s+p)

rlocus(G0)
hold on;
rlocus(G0*Gi) % mozhe da se uvidi deka GMK e skoro nepromeneto

%% proverka
clc

syms s

G0 = 20/(s*(s+10)^2)
Gi = (s+z)/(s+p)

Kv = limit(s*G0*Gi, s, 0)

e_stac = double(1/Kv)

%% 6.32 - neidealen I kompenzator
clc
clear

% so sinteza so GMK najprvo se nagoduva PREODEN rezhim, pa posle
% STACIONAREN

s = tf('s');
G0 = 1/(s*(s+5));

%rlocus(G0)
%sgrid(0.592, 50)
%ylim([-7, 7])

% od GMK se zakluchuva deka za ceta=0.592 potrebno e K = 17.8 i se dobivaat
% dominantni polovi s_1/2 = -2.5 +- j3.4

syms s K

G0 = K/(s*(s+5));

Kv = limit(s*G0, s, 0)

e_stac = 1/Kv

K_potrebno = solve(e_stac == 2)

% za dobienoto K od prethodno (za ceta):
e_stac_dobieno = 5/17.8

% dobienata vrednost na e(inf)=0.28 pri K = 17.8 e mnogu pogolema od
% potrebnata 0.02

e_stac_posakuvano = 0.02;

gain_needed = e_stac_dobieno/e_stac_posakuvano

% potrebnoto zasiluvanje vsushnost e odnosot z/p na kompenzatorot

z_p = gain_needed

z = 0.1; % proizvolno
p = z * 1/z_p;

s = tf('s');

K = 17.8
G0 = 1/(s*(s+5));
Gi = (s+z)/(s+p)

rlocus(G0);
hold on
rlocus(G0*Gi);
sgrid(0.592, 50)
ylim([-7, 7])

%% proverka
clc

syms s

G0 = K/(s*(s+5));
Gi = (s+z)/(s+p);

Kv = limit(s*G0*Gi, s, 0);

e_stac = double(1/Kv)

s = tf('s');

K = 17.8
G0 = 1/(s*(s+5));
Gi = (s+z)/(s+p)

t = 0:0.01:100;

lsim(feedback(G0, 1), t, t);
hold on;
lsim(feedback(G0*Gi, 1), t, t);