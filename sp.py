from numpy import linalg as la
from math  import e
from sympy import *

# Creates toints trigexp non-linear system
def trigexp_toint(v):
    f = [3*(v[0])**3+2*(v[1])-5+sin(v[0]-v[1])*sin(v[0]+v[1])]
    for i in range(1,len(v)-1):
        f.append(-(v[i-1])*e**(v[i-1]-v[i]) + (v[i])*(4+3*(v[i])**2)
                 + 2*(v[i+1]) + sin(v[i]-v[i+1])*sin(v[i]+v[i+1])-8)
    f.append(-(v[-2])*e**(v[-2]-v[-1])+4*v[-1]-3)

    return f

# Creates broyden tridiagonal non-linear system
def broyden_tridiagonal(v):
    f = [(3-2*v[0])*v[0]-2*v[1]+1]
    for i in range(1,len(v)-1):
        f.append((3-2*v[i])*v[i]-v[i-1]-2*v[i+1]+1)
    f.append((3-2*v[-1])*v[-1]-v[-2]+1)

    return f

# Calculates sybolicaly the jacobian of f
def jacobian (v, f):
    j = zeros(len(f),len(f))
    for l in range(0, len(f)):
        for c in range(0, len(f)):
            j[l,c] = diff(f[l], v[c])

    return j

# Calculates the solution of (j(xk)*s=f(xk)) for newtons non-linear method
def getS (j, f, x):
    fc,jc = [0]*len(f),[[0 for j in range(0, len(f))] for i in range(0, len(f))]

    for i in range(0, len(f)):
        fc[i] = float(-1*placeVars(f[i],v,x))

    for l in range(0, len(f)):
        for c in range(0, len(f)):
            jc[l][c] = float(placeVars(j[l,c],v,x))

    return la.solve(jc,fc)

# Calculates recursively f(x1,...,xn)
def placeVars (f, v, x):
    if (v != []):
        hx, *tx = x
        hv, *tv = v
        return placeVars(f.subs(hv,hx),tv,tx)

    return f

# Performs iterations to solve system by newtons method
def iterate (j,f,x,e,k):
    if (e > abs(placeVars(f[0],v,x))): return

    x += getS(j,f,x)
    print (x," ",k)
    iterate (j,f,x,e,k+1)

### MAIN ######################################################################
# variables = ['x1','x2','x3','x4','x5','x6','x7','x8','x9','x10']
variables = ['x1','x2']
v = [Symbol(variables[i]) for i in range(0,len(variables))]
f = [(v[0])**2+(v[1])**2-2,e**(v[0]-1)+(v[1])**3-2]
x = [1.5,2.0]
e = 0.0001

iterate(jacobian(v,f),f,x,e,1)
