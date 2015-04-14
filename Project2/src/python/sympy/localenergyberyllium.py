from sympy import symbols, diff, exp, sqrt, factor, Symbol, printing

import os
import numpy as np

infile = open('tmp.tex', 'w')
cfile = open('Berylliumcppcode.txt', 'w')

N, Z = 4, 4

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

psi1s = [0,0,0,0,0]
psi2s = [0,0,0,0,0]

for i in range(1,5):
    psi1s[i] = exp(-a*r[0][i])
    psi2s[i] = (1 - a*r[0][i]/2)*exp(-a*r[0][i]/2)


psiDberyllium = (psi1s[1]*psi2s[2] - psi1s[2]*psi2s[1])*\
        (psi1s[3]*psi2s[4] - psi1s[4]*psi2s[3])


def nabla(f,r,x,y,z):
    # return [0,\
    #         (diff(f,x)),\
    #         (diff(f,y)),\
    #         (diff(f,z))]

    d = diff(f,r)

    return [0,\
            (x/r)*d,\
            (y/r)*d,\
            (z/r)*d,\
            ]


g1 = nabla(psiDberyllium, r1, x1, y1, z1)[1]*(r1/x1)
g2 = nabla(psiDberyllium, r2, x2, y2, z2)[1]*(r2/x2)
g3 = nabla(psiDberyllium, r3, x3, y3, z3)[1]*(r3/x3)
g4 = nabla(psiDberyllium, r4, x4, y4, z4)[1]*(r4/x4)

G1, G2, G3, G4 = symbols('G_1 G_2 G_3 G_4')

f12 = r12*(1+b*r12)**2
f13 = r13*(1+b*r13)**2
f14 = r14*(1+b*r14)**2
f23 = r23*(1+b*r23)**2
f24 = r24*(1+b*r24)**2
f34 = r34*(1+b*r34)**2

F12, F13, F14, F23, F24, F34 =\
        symbols('F_12 F_13 F_14 F_23 F_24 F_34')

def laplace(f,r):
    return (diff(diff(f, r),r) + (2/r)*diff(f,r))

psiDberylliumX =\
    [0,\
    psiDberyllium.subs(r1, sqrt(x1*x1 + y1*y1 + z1*z1)),\
    psiDberyllium.subs(r2, sqrt(x2*x2 + y2*y2 + z2*z2)),\
    psiDberyllium.subs(r3, sqrt(x3*x3 + y3*y3 + z3*z3)),\
    psiDberyllium.subs(r4, sqrt(x4*x4 + y4*y4 + z4*z4))]

psiDnabla1 = [0,\
        nabla(psiDberyllium, r1, x1, y1, z1),\
        nabla(psiDberyllium, r2, x2, y2, z2),\
        nabla(psiDberyllium, r3, x3, y3, z3),\
        nabla(psiDberyllium, r4, x4, y4, z4)]


psiDnabla2 = [0,\
        laplace(psiDberyllium, r1),\
        laplace(psiDberyllium, r2),\
        laplace(psiDberyllium, r3),\
        laplace(psiDberyllium, r4)]

### First derivative part of jastrow factor:

psiCnabla1 = [0,[0,0,0,0],[0,0,0,0],[0,0,0,0],[0,0,0,0]]
A = 0

for dim in range(1,4):
    for k in range(1,N+1):
        for j in range(1,N+1):
            if j != k:
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
                psiCnabla1[k][dim] += (R[k][dim]-R[j][dim])\
                        *A/(r[k][j]*(1 + b*r[k][j])**2)

### Second derivative of the Jastrow factor:

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

                vecpart = (R[k][0] - R[i][0])*(R[k][0] - R[j][0])
                ijpart += vecpart/(r[k][i]*r[k][j])*\
                    (A/((1+b*r[k][i])**2))*(A/((1+b*r[k][j])**2))
    
    psiCnabla2[k] = jpart + ijpart

R1R2 = Symbol('\mathbf{R12}')

K = [0,0,0,0,0]

for i in range(1,N+1):
    tmp = 0
    for dim in range(1,4):
        tmp += psiDnabla1[i][dim]*psiCnabla1[i][dim]
    K[i] = -0.5*(psiDnabla2[i] + psiCnabla2[i] + 2*tmp)
   

def V(N, Z, r):
    V = 0
    for i in range(1,N+1):
        V -= Z/r[0][i]
        for j in range(1,N+1):
            if i<j:
                V += 1./r[i][j]
    return V

Vberyllium = V(N, Z, r)

psiSUM = sum(K) + Vberyllium

psiSUM = psiSUM.nsimplify()

psiSUM = psiSUM.subs(g1,G1).subs(g2,G2).subs(g3,G3).subs(g4,G4)

psiSUM = psiSUM\
        .subs(f12,F12).subs(f13,F13).subs(f14,F14)\
        .subs(f23,F23).subs(f24,F24).subs(f34,F34)


### C++ code printing:
C = printing.ccode(psiSUM)
cfile.write(C)
cfile.close()

### Latex printing
y1 = printing.latex(psiSUM)

x = '\documentclass[12pt,a3paper]{article}\n\
    \usepackage{amsmath}\n\
    \usepackage[landscape]{geometry}\n\
    \\begin{document}\n\
    \\tiny{$'


z = '$}\n\
    \end{document}'

infile.write(x+'\n'+y1+z)
infile.close()

os.system('pdflatex tmp.tex')
