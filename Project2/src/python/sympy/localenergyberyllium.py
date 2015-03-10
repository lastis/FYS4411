import numpy as np
from sympy import symbols, diff, exp, sqrt, factor, Symbol, printing
import os

infile = open('tmp.tex', 'w')
cfile = open('cppcodeberyllium.txt', 'w')

x1, y1, z1, x2, y2, z2, x3, y3, z3, x4, y4, z4, a, b  \
    = symbols('x_1 y_1 z_1 x_2 y_2 z_2 x_3 y_3 z_3 x_4 y_4 z_4 alpha beta')


r1 = sqrt(x1*x1 + y1*y1 + z1*z1)
r2 = sqrt(x2*x2 + y2*y2 + z2*z2)
r3 = sqrt(x3*x3 + y3*y3 + z3*z3)
r4 = sqrt(x4*x4 + y4*y4 + z4*z4)

r12 = sqrt((x1-x2)*(x1-x2) + (y1-y2)*(y1-y2) + (z1-z2)*(z1-z2))
r13 = sqrt((x1-x3)*(x1-x3) + (y1-y3)*(y1-y3) + (z1-z3)*(z1-z3))
r14 = sqrt((x1-x4)*(x1-x4) + (y1-y4)*(y1-y4) + (z1-z4)*(z1-z4))
r23 = sqrt((x2-x3)*(x2-x3) + (y2-y3)*(y2-y3) + (z2-z3)*(z2-z3))
r24 = sqrt((x2-x4)*(x2-x4) + (y2-y4)*(y2-y4) + (z2-z4)*(z2-z4))
r34 = sqrt((x3-x4)*(x3-x4) + (y3-y4)*(y3-y4) + (z3-z4)*(z3-z4))

R1 = Symbol('r_1')
R2 = Symbol('r_2')
R3 = Symbol('r_3')
R4 = Symbol('r_4')

R12 = Symbol('r_{12}')
R13 = Symbol('r_{13}')
R14 = Symbol('r_{14}')
R23 = Symbol('r_{23}')
R24 = Symbol('r_{24}')
R34 = Symbol('r_{34}')


r = [[0,r1,r2,r3,r4],[0,0,r12,r13,r14],[0,0,0,r23,r24],[0,0,0,0,r34],[]]

tmp = 0
for i in range(1,5):
    for j in range(1,5):
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
            
            tmp += A*r[i][j]/(1+b*r[i][j])

psiC = exp(tmp).nsimplify()

# .subs(r12,R12).subs(r13,R13).subs(r14,R14).subs(r23,R23)\
#         .subs(r24,R24).subs(r34,R34)\
#         .subs(r1,R1).subs(r2,R2).subs(r3,R3).subs(r4,R4)


psi1s = [0,0,0,0,0]
psi2s = [0,0,0,0,0]

for i in range(1,5):
    psi1s[i] = exp(-a*r[0][i])
    psi2s[i] = (1 - a*r[0][i]/2)*exp(-a*r[0][i]/2)


psiT = (psi1s[1]*psi2s[2] - psi1s[2]*psi2s[1])*\
        (psi1s[3]*psi2s[4] - psi1s[4]*psi2s[3])


psi = (psiT*psiC)

def laplace(f,x,y,z):
    return (diff(diff(f, x),x)+diff(diff(f, y),y)+diff(diff(f, z),z))

## Full equation
L1 = laplace(psi,x1,y1,z1).simplify().collect(psi)
L1 = (L1/psi)

L2 = laplace(psi,x2,y2,z2).simplify().collect(psi)
L2 = (L2/psi)

L3 = laplace(psi,x3,y3,z3).simplify().collect(psi)
L3 = (L3/psi)

L4 = laplace(psi,x4,y4,z4).simplify().collect(psi)
L4 = (L4/psi)

EL = -L1/2 - L2/2 -L3/2 -L4/2 - 4/r1 - 4/r2 - 4/r3 - 4/r4 \
        + 1/r12 + 1/r13 + 1/r14 + 1/r23 + 1/r24 + 1/r34


x = '\documentclass[12pt,a3paper]{article}\n\
    \usepackage{amsmath}\n\
    \usepackage[landscape]{geometry}\n\
    \\begin{document}\n\
    \\tiny{$'


EL = EL.subs(r12,R12).subs(r13,R13).subs(r14,R14).subs(r23,R23)\
        .subs(r24,R24).subs(r34,R34)\
        .subs(r1,R1).subs(r2,R2).subs(r3,R3).subs(r4,R4).simplify()\
        .subs(r12**2,R12**2).subs(r13**2,R13**2)\
        .subs(r14**2,R14**2).subs(r23**2,R23**2)\
        .subs(r24**2,R24**2).subs(r34**2,R34**2)\
        .subs(r1**2,R1**2).subs(r2**2,R2**2)\
        .subs(r3**2,R3**2).subs(r4**2,R4**2)


y = printing.latex(EL)

c = printing.ccode(EL)

z = '$}\n\
    \end{document}'

infile.write(x+'\n'+y+z)
infile.close()

cfile.write(c)
cfile.close()

os.system('pdflatex tmp.tex')

