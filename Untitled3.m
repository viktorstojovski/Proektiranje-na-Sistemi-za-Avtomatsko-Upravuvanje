syms K T
diffvar = {'K','T'};
numden{1} = [(T-K) (1-K)];
numden{2} = [T 1 K];
JdJ_funkcija(numden,diffvar);