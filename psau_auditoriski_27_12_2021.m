%% Prva (proektiranje na sostojben upravuvach)
clc
clear

clc
clear

Mp = 20.8;
C = (-log(Mp/100))/sqrt(pi^2 + (log(Mp/100))^2)
wn = 2.23

s12 = -C*wn + sqrt(-1)*wn*sqrt(1-C^2)


A = [-5 1 0; 0 -2 1; 0 0 -1];
B = [0 0 1]';
C = [-1 1 0];
D = 0

K = [0 5 -1];

A_c = A - B*K

[num, den] = ss2tf(A, B, C, D)
[num_c, den_c] = ss2tf(A_c, B, C, D)

G_nc = tf(num, den)
G_c = tf(num_c, den_c)

pole(G_c)

step(G_nc)
hold on;
step(G_c)
legend('Nekopmenziran', 'Kompenziran')

%% Vtora (sinteza preku upravliva kanonska forma)
clc
clear

Mp = 20.8;
C = (-log(Mp/100))/sqrt(pi^2 + (log(Mp/100))^2)
wn = 2.23

s12 = -C*wn + sqrt(-1)*wn*sqrt(1-C^2)

A = [-5 1 0; 0 -2 1; 0 0 -1]
B = [0 0 1]'
C = [-1 1 0]
D = [0]

sys = ss(A,B,C,D)

[num, den] = ss2tf(A, B, C, D)

G_nc = tf(num, den)

syms s kv1 kv2 kv3

Av = [0 1 0; 0 0 1; -10 -17 -8];
Bv = [0 0 1]';

Kv = [kv1 kv2 kv3];

as = collect((s+1-1j*2)*(s+1+1j*2)*(s+5));

akv = collect(det(s*eye(3) - (Av - Bv*Kv)))

sol = solve(coeffs(akv, s) == coeffs(as, s), [kv1 kv2 kv3]);

Kv = [sol.kv1, sol.kv2, sol.kv3]

%% naogjanje na P
clc

Qcw = ctrb(A, B)

fprintf('Proverka na upravlivsot na originalniot sistem:\n')
if rank(Qcw) == length(A)
   fprintf('Rangot na Qcw e %d, dimenziite na A se %dx%d => originalniot sistem e kompletno upravliv\n', rank(Qcw), length(A), length(A))
else
   fprintf('Rangot na Qcw e %d, dimenziite na A se %dx%d => originalniot sistem NE e kompletno upravliv\n', rank(Qcw), length(A), length(A))
end

Qcv = ctrb(Av, Bv)

fprintf('Proverka na upravlivsot na transformiraniot sistem:\n')
if rank(Qcv) == length(Av)
   fprintf('Rangot na Qcv e %d, dimenziite na A se %dx%d => transformiraniot sistem e kompletno upravliv\n', rank(Qcw), length(A), length(A))
else
   fprintf('Rangot na Qcv e %d, dimenziite na A se %dx%d => transformiraniot sistem NE e kompletno upravliv\n', rank(Qcw), length(A), length(A))
end

P = Qcw * inv(Qcv)

Kw = Kv * inv(P) % se poklopuva so dobienoto K od zadacha 1!

%% Treta (proektiranje na observer preku nabljudliva kanonska forma)
clc
clear

syms s

as = collect((s+1)^2 * (s+2))
aso = collect((s+10)^2 * (s+20))

A = [0 1 0; 0 0 1; -4 -6 -4];
B = [0 0 1]';
C = [1 0 1];
D = 0;

syms k1 k2 k3

K = [k1 k2 k3]

ak = collect(det(s*eye(3) - (A-B*K)))

sol = solve(coeffs(ak, s) == coeffs(as, s), [k1 k2 k3]);

K = [sol.k1, sol.k2, sol.k3]

Av = A.';
Bv = C.';
Cv = B.';
Dv = D;

syms Lv1 Lv2 Lv3

Lv = [Lv1 Lv2 Lv3].';

alvo = collect(det(s*eye(3) - (Av-Lv*Cv)))

sol = solve(coeffs(alvo, s) == coeffs(aso, s), [Lv1 Lv2 Lv3]);

Lv = [sol.Lv1, sol.Lv2, sol.Lv3].'

%% naogjanje na P
clc

Qow = obsv(A, C)

fprintf('Proverka na nabljudlivost na originalniot sistem:\n')
if rank(Qow) == length(A)
   fprintf('Rangot na Qcw e %d, dimenziite na A se %dx%d => originalniot sistem e kompletno nabljudliv\n', rank(Qow), length(A), length(A))
else
   fprintf('Rangot na Qcw e %d, dimenziite na A se %dx%d => originalniot sistem NE e kompletno nabljudliv\n', rank(Qow), length(A), length(A))
end

Qov = obsv(Av, Cv)

fprintf('Proverka na nabljudlivost na transformiraniot sistem:\n')
if rank(Qov) == length(Av)
   fprintf('Rangot na Qcv e %d, dimenziite na A se %dx%d => transformiraniot sistem e kompletno nabljudliv\n', rank(Qow), length(A), length(A))
else
   fprintf('Rangot na Qcv e %d, dimenziite na A se %dx%d => transformiraniot sistem NE e kompletno nabljudliv\n', rank(Qow), length(A), length(A))
end

P = inv(Qow) * Qov

Lw = vpa(P * Lv, 5)

polovi_observer = vpa(eig(A - Lw*C), 5)
polovi_sistemot = vpa(eig(A - B*K), 5)