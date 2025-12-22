import numpy as np
P=700000
T=313.15
R=8.314
Pc=1.1353*10**7
Tc=405.4
ω=0.256
a=(0.45724*R**2*T**2)/Pc
b=(0.07780*R*T)/Pc
#équation
#P*Vm**3+(b*P-R*T)*Vm**2+Vm*(a-3*P*b**2-2*R*T*b)-(R*T*b**2-b*a+P*b**3)=0
#coefficients
A=P
B=(b*P-R*T)
C=(-2*R*T*b+a-3*P*b**2)
D=(R*T*b**2-b*a+P*b**3)
racine=np.roots([A,B,C,D])
print(racine)