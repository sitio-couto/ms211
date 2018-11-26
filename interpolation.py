from sympy import *
from numpy import matrix, zeros, linalg, inf

def print_diff_table(orders):
    for i in range(len(orders)):
        print("Order", i)
        for x in orders[i]:
            print("%.2f" % x)
        print("\n")

def diff_table(x, y):
    dim = 2*len(y)-1
    table = [[inf for i in range(dim)] for y in range(len(y))]
    xt = [inf for i in range(dim)]

    for i in range(0,dim,2):
        table[0][i] = y[int(i/2)]
        xt[i] = x[int(i/2)]

    for ord in range(1, len(y)):
        prev = table[ord-1]
        for i in range(ord, dim-ord, 2):
            table[ord][i] = (prev[i+1]-prev[i-1])/(xt[i+ord]-xt[i-ord])

    print("\nDiff table")
    for i in range(len(y)): print("%-10s" %("Order "+str(i)), end="")
    for j in range(dim):
        print("")
        for i in range(len(y)):
            if table[i][j] == inf : print("%-10s" %"", end="")
            elif table[i][j] >= 0 : print(" %-9.2f" %table[i][j], end="")
            else : print("%-10.2f" %table[i][j], end="")

    print("")
    return

def vandermonde(x, y):
    dim = len(x)
    a = matrix(zeros((dim, dim)))
    for i in range(dim):
        for j in range(dim):
            a[i,j] = x[i]**j

    return linalg.solve(a, y)

def lagrange(var, x, y):
    pn = 0

    for i in range(len(y)):
        lx = 1
        for j in range(len(y)):
            if j != i: lx *= (var - x[j])/(x[i] - x[j])
        pn += y[i]*lx

    return pn

def bissect(v, f, i, e):
    x = (i[0]+i[1])/2
    y = f.subs(v, x)

    if y < 0 : i[0] = x
    else: i[1] = x

    if (i[1]-i[0]) < e : return x
    else: return bissect(v, f, i, e)

def min_squares(var, g, x, y):
    dim = len(g)
    a = matrix(zeros((dim, dim)))
    b = [0]*dim

    for i in range(dim):
        for j in range(len(y)):
            b[i] += float(g[i].subs(var, x[j])*y[j])

    for i in range(dim):
        for j in range(dim):
            for k in x :
                a[i,j] += g[i].subs(var, k)*g[j].subs(var, k)

    return linalg.solve(a, b)

def euler_improved(x, y, x0, y0, f, h, end):
    if x0 >= end : return []
    aux = f.subs([(x,x0),(y,y0)])
    y1 = y0 + (h/2)*(aux + f.subs([(x,x0+h),(y,y0+h*aux)]))
    return [y1]+euler_improved(x, y, x0+h, y1, f, h, end)

#### ENTRADAS E SAIDAS #########################################################

# Exercício 1
var1 = Symbol("x")
x1 = [-4.2,-1.4,0.8,2]
y1 = [1.7,2.5,3,2]
print("\nSAÍDA EXERCÍCIO 1")
pn = lagrange(var1, x1[:3], y1[:3])
print(pn.subs(var1, -3))
diff_table([-4.2,-1.4,0.8,2], [1.7,2.5,3,2])

#Exercício 2
h2 = Symbol("h")
a2 = Symbol("Xi")
b2 = Symbol("Xi+1")
c2 = Symbol("Xi+2")
x2 = Symbol("x")
f2 = x2**3 + ln(x2)
sp = (h2/3)*(f2.subs(x2,a2) + 4*f2.subs(x2,b2) + f2.subs(x2,c2))
print("\nSAÍDA EXERCÍCIO 2")
step1 = [(h2,0.25),(a2,1),(b2,1.25),(c2,1.5)]
step2 = [(h2,0.25),(a2,1.5),(b2,1.75),(c2,2)]
result = sp.subs(step1) + sp.subs(step2)
print(result.evalf())
# Exercício 3
var3 = Symbol("x")
x3 = [-2, -1, 0, 1, 2]
y3 = [x**-2 for x in [1, 2, 3, 2, 1]]
g3 = [var3**2, sympify(1)]
print("\nSAÍDA EXERCÍCIO 3")
print(min_squares(var3, g3, x3, y3))

# Exercício 4
x = Symbol("X")
a = Symbol("Yi-1")
b = Symbol("Yi")
c = Symbol("Yi+1")
h = Symbol("h")
dy1 = (c-a)/(2*h)
dy2 = (c - 2*b + a)/h**2
deq = dy2 + x*dy1 + b - 5*x
print("\nSAIDA EXERCÍCIO 4")
print(deq.subs([(x,0.25),(a,1),(h,0.25)]))
print(deq.subs([(x,0.5),(h,0.25)]))
print(deq.subs([(x,0.75),(c,0),(h,0.25)]))

# Entrada exercício 5
var5 = Symbol("t")
x5 = [0, 1, 2, 3, 4]
y5 = [200, 400, 650, 850, 950]
g5 = [sympify(1), var5]
print("\nSAÍDA EXERCÍCIO 5 (A)")
print(min_squares(var5, g5, x5, [log(1000/x - 1) for x in y5]))
print("\nSAÍDA EXERCÍCIO 5 (C)")
print(bissect(var5, lagrange(var5, x5[:3], y5[:3])-480, [1, 1.5], 0.001))
print("\nSAÍDA EXERCÍCIO 5 (D)")
var5_yt = Symbol("y")
dy = -var5_yt*(var5_yt/1000 - 1)
print(euler_improved(Symbol("any"), var5_yt, 0, 200, dy, 0.5, 1))
