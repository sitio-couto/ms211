from sympy import *
from numpy import matrix, zeros, linalg

x1 = [1.7,2,2.5,3]
y1 = [-4.2,2,-1.4,0.8]

var3 = Symbol("x")
x3 = [-2, -1, 0, 1, 2]
y3 = [1, 1/4, 1/9, 1/4, 1]
g3 = [var3**2, sympify(1)]


def print_diff_table(orders):
    for i in range(len(orders)):
        print("Order", i)
        for x in orders[i]:
            print("%.2f" % x)
        print("\n")

def diff_table(x, y):
    orders = [y.copy()]

    for ord in range(1,len(x)):
        new_ord = []
        for i in range(len(x)-ord):
            last_ord = orders[-1]
            new_ord.append((last_ord[i+1]-last_ord[i])/(x[i+ord]-x[i]))
        orders.append(new_ord.copy())

    return orders

def vandermonde(x, y):
    dim = len(x)
    a = matrix(zeros((dim, dim)))
    for i in range(dim):
        for j in range(dim):
            a[i,j] = x[i]**j

    return linalg.solve(a, y)

def lagrange(x, y):
    pn = 0
    var = Symbol("x")

    for i in range(len(y)):
        lx = 1
        for j in range(len(y)):
            if j != i: lx *= (var - x[j])/(x[i] - x[j])
        pn += y[i]*lx

    return pn

def bissect(f, i, e):
    x = (i[0]+i[1])/2
    y = f.subs(Symbol("x"), x)

    if y < 0 : i[0] = x
    else: i[1] = x

    if (i[1]-i[0]) < e : return x
    else: return bissect(f, i, e)

def min_squares(var, g, x, y):
    dim = len(g)
    a = matrix(zeros((dim, dim)))
    b = [float(y[i]*g[i].subs(var, x[i])) for i in range(dim)]

    for i in range(dim):
        for j in range(dim):
            for k in x :
                a[i,j] += g[i].subs(var, k)*g[j].subs(var, k)

    return linalg.solve(a, b)

def finite_diferences()

pn = lagrange(x1[:3], y1[:3])
print(bissect(pn + 3, [1.7, 2], 0.01))

print(min_squares(var3, g3, x3, y3))
# print_diff_table(diff_table(x, y))
