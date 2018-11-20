from sympy import *
from numpy import matrix, zeros, linalg, log

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
    b = [0]*dim

    for i in range(dim):
        for j in range(len(y)):
            b[i] += float(g[i].subs(var, x[j])*y[j])

    for i in range(dim):
        for j in range(dim):
            for k in x :
                a[i,j] += g[i].subs(var, k)*g[j].subs(var, k)

    return linalg.solve(a, b)

#### ENTRADAS E SAIDAS #########################################################

# Entrada exercício 1
x1 = [1.7,2,2.5,3]
y1 = [-4.2,2,-1.4,0.8]
# print("SAÍDA EXERCÍCIO 1")
# pn = lagrange(x1[:3], y1[:3])
# print(bissect(pn + 3, [1.7, 2], 0.01))

# Entrada exercício 3
var3 = Symbol("x")
x3 = [-2, -1, 0, 1, 2]
y3 = [x**-2 for x in [1, 2, 3, 2, 1]]
g3 = [var3**2, sympify(1)]
print("SAÍDA EXERCÍCIO 3")
print(min_squares(var3, g3, x3, y3))

# Entrada exercício 5
var5 = Symbol("t")
x5 = [0, 1, 2, 3, 4]
y5 = [log(1000/x - 1) for x in [200, 400, 650, 850, 950]]
g5 = [sympify(1), var5]
print("SAÍDA EXERCÍCIO 5")
print(min_squares(var5, g5, x5, y5))

# # Teste de corretude (baseado no exemplo de quadrados minimos nao linear do livro)
# varex = Symbol("x")
# xex = [-1, -0.7, -0.4, -0.1, 0.2, 0.5, 0.8, 1]
# yex = [log(x) for x in [36.547, 17.264, 8.155, 3.852, 1.820, 0.860, 0.406, 0.246]]
# gex = [sympify(1), varex]
# print("EXEMPLO:")
# print(min_squares(varex, gex, xex, yex))
x = Symbol("X")
a = Symbol("Yi-1")
b = Symbol("Yi")
c = Symbol("Yi+1")
h = Symbol("h")
dy1 = (c-a)/(2*h)
dy2 = (c - 2*b + a)/h**2
deq = dy2 + x*dy1 + b - 5*x
