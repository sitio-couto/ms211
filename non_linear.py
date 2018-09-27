import math as m
import numpy as np
import symengine as se

x = [-1,-1,-1,-1,-1,-1,-1,-1,-1,-1]
e = 0.0001

def fa(x,y,z): return (3-2*x)*x-2*y+1
def fb(x,y,z): return (3-2*x)*x-z-2*y+1
def fc(x,y,z): return (3-2*x)*x-z+1
def fal (x): return -2*x + 3 - 2*x


def f1(x,y): return (3-2*x)*x-2*y+1
def f2(x,y): return (3-2*x)*x-2*y+1

def f1x (x): return -20*x
def f1y (y): return 10
def f2x (x): return -1
def f2y (y): return 0

def iteratesp (x,e,k):
    f = [-fa (x[0],x[1],0)]
    for i in range(1,9):
        f.append(-fb (x[i], x[i+1], x[i-1]))
    f.append(-fc (x[-1], 0, x[-2]))
    if (abs(f[0]) < e): return x.copy(), k

    j =[[0 for x in range(10)] for y in range(10)]
    j[0] = [fal(x[0]),-2, 0, 0,0,0,0,0,0,0]
    for i in range(1,9):
        j[i][i-1], j[i][i], j[i][i+1] = fal(x[i]), -2, -1
    j[-1] = [ 0,0,0,0,0,0,0, fal(x[-1]), 0, -1]

    s = np.linalg.solve(j,f)
    x = np.add(x,s)

    print (x)
    return iteratesp (x,e,k+1)

def iterate (x,e,k):
    f = [-f1(x[0], x[1]), -f2(x[0], x[1])]
    if (abs(f[0]) < e): return x.copy(), k
    j = [[f1x(x[0]),f1y(x[1])],[f2x(x[0]),f2y(x[1])]]
    s = np.linalg.solve(j,f)
    x = np.add(x,s)
    print (x)
    return iterate (x,e,k+1)

a, b = iteratesp(x, e, 0)
# a, b = iterate (x,e,0)
# print (a," ",b)

vars = 'x y z'
func = ['(3-2*x)*x-2*y+1','(3-2*x)*x-z-2*y+1','(3-2*x)*x-z+1']
v,f = se.symbols(vars),[]
for i in range(0, len(func)): f.append(se.sympify(func[i]))

def jacob (v, f):
    jf = se.zeros(len(f),len(v))
    for i in range(0, len(f)):
        for j in range(0, len(v)):
            jf[i,j] = se.diff(f[i], v[j])

    return jf

# print (jacob(v,f))
# def iterate ():
