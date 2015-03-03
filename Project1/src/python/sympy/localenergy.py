from sympy import symbols, diff, exp, sqrt, factor, Symbol, printing
import os

infile = open('tmp.tex', 'w')

x1, y1, z1, x2, y2, z2, a, b  = symbols('x_1 y_1 z_1 x_2 y_2 z_2 a b')

r1 = sqrt(x1*x1 + y1*y1 + z1*z1)
r2 = sqrt(x2*x2 + y2*y2 + z2*z2)

r12 = sqrt((x1-x2)*(x1-x2)+ (y1-y2)*(y1-y2)+ (z1-z2)*(z1-z2))

rdot = x1*x2 + y1*y2 + z1*z2

psiT1 = exp(-a*(r1+r2))
f = b*r12+1
psiT2 = exp(-a*(r1+r2) + r12/(2*f))

x,y,z = x1-x2, y1-y2, z1-z2
X,Y,Z = symbols('x y z')

r12 = sqrt(x**2 + y**2 + z**2)

R1 = Symbol('r_1')
R2 = Symbol('r_2')
R12 = Symbol('r_{12}')
Rdot = Symbol('\mathbf{r}')
F = Symbol('f')

def laplace(f,x,y,z):
    return (diff(diff(f, x),x)+diff(diff(f, y),y)+diff(diff(f, z),z))

def dd(f,x):
    return diff(diff(f,x),x)

## Double derivative for x_1 - WORKS! 
'''
B = dd(psiT2,x1).collect(psiT2)
B = (B/psiT2)
B = B.subs(r1,R1).subs(r2,R2).subs(r1**2,R1**2).subs(r2**2,R2**2).subs(r12,R12)
B = B.simplify()
print printing.latex(B)
'''

## Laplace for x_1 and x_2, summed and *(-1/2) - WORKS!
'''
L1 = dd(psiT2,x1).collect(psiT2)
L1 = (L1/psiT2)
L1 = L1.subs(r1,R1).subs(r2,R2).subs(r1**2,R1**2)\
        .subs(r2**2,R2**2).subs(r12,R12)

L2 = dd(psiT2,x2).collect(psiT2)
L2 = (L2/psiT2)
L2 = L2.subs(r1,R1).subs(r2,R2).subs(r1**2,R1**2)\
        .subs(r2**2,R2**2).subs(r12,R12)

L = (-L1/2-L2/2).simplify()

print printing.latex(L)
'''

## Full equation
L1 = laplace(psiT2,x1,y1,z1).collect(psiT2)
L1 = (L1/psiT2)

L2 = laplace(psiT2,x2,y2,z2).collect(psiT2)
L2 = (L2/psiT2)

EL2 = -L1/2-L2/2-2/r1-2/r2+1/r12

## Kinda works:
'''
EL2 = EL2.subs(r1,R1).subs(r2,R2).subs(r12,R12).simplify()\
        .subs(r1**2,R1**2).subs(r2**2,R2**2).subs(r12**2,R12**2)\
        .subs(f.subs(r12,R12),F)

EL2 = EL2.expand().simplify()


x = '\documentclass[12pt,a3paper]{article}\n\
    \usepackage{amsmath}\n\
    \usepackage[landscape]{geometry}\n\
    \\begin{document}\n\
    \\tiny{$'

el1 = a/r1 + a/r2 - 2/r1 - 2/r2 + 1/r12 - a**2

el1 = el1.subs(r1,R1).subs(r2,R2).subs(r12,R12)

EL1 = Symbol('E_{L1}')

EL2 = EL2.subs(r1**2,R1**2).subs(r2**2,R2**2)\
        .subs(r12**2,R12**2).subs(rdot, Rdot).subs(2*rdot, 2*Rdot).expand()

EL2 = EL2.simplify().subs(el1, EL1)

y = printing.latex(EL2.subs(F,f.subs(r12,R12)))
'''

x = '\documentclass[12pt,a3paper]{article}\n\
    \usepackage{amsmath}\n\
    \usepackage[landscape]{geometry}\n\
    \\begin{document}\n\
    \\tiny{$'


EL2 = EL2.subs(r1,R1).subs(r2,R2).subs(r12,R12)\
        .expand().simplify().subs(r12**2,R12**2)\
        .subs(r1**2,R1**2).subs(r2**2,R2**2).subs(rdot,Rdot)\
        .subs(rdot**2,Rdot**2)

y = printing.latex(EL2)

z = '$}\n\
    \end{document}'

infile.write(x+'\n'+y+z)
infile.close()

os.system('pdflatex tmp.tex')
print '\n', 1/((2*f**2).subs(r12,R12))
