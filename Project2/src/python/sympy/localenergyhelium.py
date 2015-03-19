from sympy import symbols, diff, exp, sqrt, factor, Symbol, printing
import os
import numpy as np

infile = open('tmp.tex', 'w')

N = 2

a, b = symbols('alpha beta')
r1, r2, r3, r4, r12, r13, r14, r23, r24, r34 \
        = symbols('r_1 r_2 r_3 r_4 r_12 r_13 r_14 r_23 r_24 r_34')

R1, R2, R3, R4 = symbols('R_1, R_2, R_3, R_4')

R = [0, R1, R2, R3, R4]

r = [[0,r1,r2,r3,r4],[0,0,r12,r13,r14],\
        [0,0,0,r23,r24],[0,0,0,0,r34],[0,0,0,0,0]]

A = 1./2

for i in range(1,N+1):
    for j in range(1,N+1):
        if i < j:
            if i == 1 and j == 2:
                A = 1./4
            elif i == 1 and j == 3:
                A = 1./2
            elif i == 1 and j == 4:
                A = 1./2
            elif i == 2 and j == 3:
                A = 1./2
            elif i == 2 and j == 4:
                A = 1./2
            elif i == 3 and j == 4:
                A = 1./4

psiCnabla = [0,0,0,0,0]

vecpart = 0

for k in range(1,N+1):
    jpart = 0
    ijpart = 0
    for j in range(1,N+1):
        if k < j:
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
            if k < j and k < i:
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
                A = 1./2
                vecpart = (R[k] - R[i])*(R[k] - R[j])
                ijpart += vecpart/(r[k][i]*r[k][j])*\
                    (A/((1+b*r[k][i])**2))*(A/((1+b*r[k][j])**2))
    
    psiCnabla[k] = jpart + ijpart

psiSUM = sum(psiCnabla).nsimplify()

R1R2 = Symbol('\mathbf{R_12}')

psiSUM = psiSUM.subs((R1-R2)**2, r1**2 + r2**2 -2*R1R2)

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
