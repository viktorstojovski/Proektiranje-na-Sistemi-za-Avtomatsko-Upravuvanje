function [J]=JdJ_funkcija(NumDen,diffvar)

% num - vector of real and symbolical values
% den - vector of real and symbolical values
% diffvar - cell
% type - J1, J2, J3
cnt = 0;
Jv=0;
for i = 2:2:length(NumDen)
    num = NumDen{i-1};
    den = NumDen{i};

if(length(den)==3)

a2 = den(1);
a1 = den(2);
a0 = den(3);
b1 = num(1);
b0 = num(2);

J = (a0*b1*b1+a2*b0*b0)/(2*a0*a1*a2);
J = simplify(J);
eval(['Jout' num2str(cnt+1) '=' 'J'])
end

if(length(den)==4)
    
a3 = den(1);
a2 = den(2);
a1 = den(3);
a0 = den(4);
b2 = num(1);
b1 = num(2);
b0 = num(3);

J = (a0*a1*b2*b2+a0*a3*(b1*b1-2*b0*b2)+a2*a3*b0*b0)/(2*a0*a3*(a1*a2-a0*a3));
J = simplify(J);
eval(['Jout' num2str(cnt+1) '=' 'J'])

end

if(length(den)==5)

a4 = den(1);
a3 = den(2);
a2 = den(3);
a1 = den(4);
a0 = den(5);
b3 = num(1);
b2 = num(2);
b1 = num(3);
b0 = num(4);

J1 = (b3*b3*(a0*a1*a2-a0*a0*a3)+a0*a1*a4*(b2*b2-2*b1*b3)+a0*a3*a4*(b1*b1-2*b0*b2))/(2*a0*a4*(a1*a2*a3-a0*a3*a3-a1*a1*a4));
J2 = (b0*b0*(a2*a3*a4-a1*a4*a4))/(2*a0*a4*(a1*a2*a3-a0*a3*a3-a1*a1*a4));
J = J1+J2;
J = simplify(J);
eval(['Jout' num2str(cnt+1) '=' 'J'])   
end
Jv = Jv+J;
cnt = cnt+1;
end
if(cnt~=1)
eval(['Jtotal' '=' 'simplify(collect(Jv))']) 
end
    for i = 1:length(diffvar)
  eval(['diffJ' diffvar{i} '=' 'simplify(collect(diff(Jv,diffvar{i}),diffvar{i}))'])
    end
    
    
    
end

