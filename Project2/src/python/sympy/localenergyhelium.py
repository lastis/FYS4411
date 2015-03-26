from sympy import symbols, diff, exp, sqrt, factor, Symbol, printing, collect
import os
import numpy as np

infile = open('tmp.tex', 'w')

N = 2

a, b = symbols('alpha beta')

x1,y1,z1,x2,y2,z2,x3,y3,z3,x4,y4,z4 \
    = symbols('x_1 y_1 z_1 x_2 y_2 z_2 x_3 y_3 z_3 x_4 y_4 z_4')

r1, r2, r3, r4, r12, r13, r14, r23, r24, r34 \
        = symbols('r_1 r_2 r_3 r_4 r_12 r_13 r_14 r_23 r_24 r_34')

R1, R2, R3, R4 = symbols('R_1, R_2, R_3, R_4')

R = [0, [R1,x1,y1,z1], [R2,x2,y2,z2], [R3,x3,y3,z3], [R4,x4,y4,z4]]

r = [[0,r1,r2,r3,r4],[0,0,r12,r13,r14],\
        [0,r12,0,r23,r24],[0,r13,r23,0,r34],[0,r14,r24,r34,0]]


### Slater determinant part and its derivatives:

psiDhelium = exp(-a*(r1+r2)) # Helium psiD

psi1s = [0,0,0,0,0]
psi2s = [0,0,0,0,0]

for i in range(1,5):
    psi1s[i] = exp(-a*r[0][i])
    psi2s[i] = (1 - a*r[0][i]/2)*exp(-a*r[0][i]/2)


psiDberyllium = (psi1s[1]*psi2s[2] - psi1s[2]*psi2s[1])*\
        (psi1s[3]*psi2s[4] - psi1s[4]*psi2s[3])


def nabla(f,x,y,z):
    return [0, (diff(f,x))/f, (diff(f,y))/f, (diff(f,z))/f]

def laplace(f,r):
    return (diff(diff(f, r),r) + (2/r)*diff(f,r)).collect(f)/f

# psiDnabla1 = [0, nabla(psiDhelium,r1), nabla(psiDhelium,r2)]


psiDheliumX =\
    [0,psiDhelium.subs(r1,sqrt(x1*x1 + y1*y1 + z1*z1)),\
    psiDhelium.subs(r2,sqrt(x2*x2 + y2*y2 + z2*z2))]

psiDnabla1 = [0,nabla(psiDheliumX[1],x1,y1,z1),nabla(psiDheliumX[2],x2,y2,z2)]
psiDnabla2 = [0, laplace(psiDhelium,r1), laplace(psiDhelium,r2)]

### First derivative part of jastrow factor:

psiCnabla1 = [0,[0,0,0,0],[0,0,0,0],[0,0,0,0],[0,0,0,0]]
A = 0

for dim in range(1,4):
    for k in range(1,N+1):
        for j in range(1,N+1):
            if j != k:
                # if k == 1 and j == 2:
                #     A = 1./4
                # elif k == 1 and j == 3:
                #     A = 1./2
                # elif k == 1 and j == 4:
                #     A = 1./2
                # elif k == 2 and j == 3:
                #     A = 1./2
                # elif k == 2 and j == 4:
                #     A = 1./2
                # elif k == 3 and j == 4:
                #     A = 1./4
                A = 1./2 ## FOR HELIUM / BERYLLIUM
                psiCnabla1[k][dim] += (R[k][dim]-R[j][dim])\
                        *A/(r[k][j]*(1 + b*r[k][j])**2)

### Double derivative part of Jastrow factor:
psiCnabla2 = [0,0,0,0,0]

vecpart = 0

for k in range(1,N+1):
    jpart = 0
    ijpart = 0
    for j in range(1,N+1):
        if k != j:
            if k == 1 and j == 2:
                A = 1./4
            elif k == 1 and j == 3:
                A = 1./2
            elif k == 1 and j == 4:
                A = 1./2
            elif k == 2 and j == 3:
                A = 1./2
            elif k == 2 and j == 4:
                A = 1./2
            elif k == 3 and j == 4:
                A = 1./4
            
            A = 1./2 ## FOR HELIUM / BERYLLIUM
            jpart += 2*A/(r[k][j]*(1+b*r[k][j])**2) -\
                2*A*b/((1+b*r[k][j])**3)
    
    for i in range(1,N+1):
        for j in range(1, N+1):
            if k != j and k != i:
                if k == 1 and j == 2:
                    A = 1./4
                elif k == 1 and j == 3:
                    A = 1./2
                elif k == 1 and j == 4:
                    A = 1./2
                elif k == 2 and j == 3:
                    A = 1./2
                elif k == 2 and j == 4:
                    A = 1./2
                elif k == 3 and j == 4:
                    A = 1./4
                A = 1./2 ## FOR HELIUM / BERYLLIUM
                vecpart = (R[k][0] - R[i][0])*(R[k][0] - R[j][0])
                ijpart += vecpart/(r[k][i]*r[k][j])*\
                    (A/((1+b*r[k][i])**2))*(A/((1+b*r[k][j])**2))
    
    psiCnabla2[k] = jpart + ijpart



R1R2 = Symbol('\mathbf{R_12}')

K = [0,0,0,0,0]

for i in range(1,N+1):
    tmp = 0
    for dim in range(1,4):
        tmp += psiDnabla1[i][dim]*psiCnabla1[i][dim]
    K[i] = -0.5*(psiDnabla2[i] + psiCnabla2[i] + 2*tmp)
   
el1 = a/r1 + a/r2 - 2/r1 - 2/r2 + 1/r12 - a**2

EL1 = Symbol('E_{L1}')

Vhelium = -2*(1/r1 + 1/r2) + 1/r12

psiSUM = sum(K) + Vhelium
# psiSUM = psiSUM.subs((R1-R2)**2, r1**2 + r2**2 -2*R1R2)
psiSUM = psiSUM.nsimplify().subs(el1,EL1)
psiSUM = psiSUM.subs(sqrt(x1**2 + y1**2 + z1**2),r1)\
        .subs(sqrt(x2**2 + y2**2 + z2**2),r2)

psiSUM = psiSUM - EL1

# psiSUM = psiSUM.simplify()


psiSUM = psiSUM.subs((R1-R2)**2, r12**2).subs((-R1+R2)**2, r12**2)

psiSUM = collect(psiSUM,\
    (a/(r12*r2*(b*r12 + 1)**2), a/(r12*r1*(b*r12 + 1)**2)))


psiSUM = psiSUM\
    .subs(x2*(-x1/2+x2/2) + y2*(-y1/2+y2/2) + z2*(-z1/2+z2/2), -R1R2 + r2**2)

psiSUM = psiSUM\
        .subs(x1*(x1/2-x2/2) + y1*(y1/2-y2/2) + z1*(z1/2-z2/2), -R1R2 + r1**2)

y = printing.latex(psiSUM)

x = '\documentclass[12pt,a3paper]{article}\n\
    \usepackage{amsmath}\n\
    \usepackage[landscape]{geometry}\n\
    \\begin{document}\n\
    \\tiny{$'


z = '$}\n\
    \end{document}'

infile.write(x+'\n'+y+z)
infile.close()

os.system('pdflatex tmp.tex')
