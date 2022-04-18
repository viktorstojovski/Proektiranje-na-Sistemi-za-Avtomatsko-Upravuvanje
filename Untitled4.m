z=tf('z')
s=tf('s')

P_s=1/(s*(s+1)*(s+2))
P_z=c2d(P_s,0.1,'zoh')
zpk(P_z)





