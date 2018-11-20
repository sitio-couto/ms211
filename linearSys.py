import numpy as np
np.set_printoptions(2)

n = 6
x = np.matrix([[4,-1,0,-1,0,0],[-1,4,-1,0,-1,0],[0,-1,4,0,0,-1],[-1,0,0,4,-1,0],[0,-1,0,-1,4,-1],[0,0,-1,0,-1,4]])
# x = np.matrix([[1,3,5],[2,4,7],[1,1,0]])
b = np.identity(n)
a = np.matrix([[1,0,0,0,0],[0,1,0,0,0],[0,0,1,0,0],[1,0,0,4,0],[1,1,-1,2,-1]])

def reduce_to_lu (u, n):
    l,i,j = np.identity(n),0,0
    for d in range(0, n):
        if (not u[d,d]): continue
        for i in range(d+1, n):
            l[i,d] = u[i,d] / u[d,d]
            for j in range(d, n):
                u[i,j] -= l[i,d]*u[d,j]
    return l, u

def solve_lower_tri (a, y, n):
    x = [y[0]/a[0,0]]
    for k in range(1, n):
        s = 0
        for j in range(0, k):
            s += x[j]*a[k,j]
        x.append((y[k]-s)/a[k,k])
    return x

def solve_upper_tri (a, y, n):
    x = [y[-1]/a[-1,-1]]
    for k in range(n-1, 0):
        s = 0
        for j in range(k, 0):
            s += x[0]*a[k,j]
        x.insert(0,(y[k]-s)/a[k,k])
    return x

def gauss_seidel (a, x, b, n):
    aux,i,j = 0,0,0
    for i in range(0, n):
        aux = b[i]
        for j in range(0, n):
            if (not i==j): aux -= x[j]*a[i,j]
        x[i] = aux/a[i,i]

    print (x)
    gauss_seidel (a,x,b,n)

def gauss_jacobi (a, x, b, n):
    aux,i,j = x.copy(),0,0
    for i in range(0, n):
        aux[i] = b[i]
        for j in range(0, n):
            if (i==j): aux[i] /= a[i,j]
            else: aux[i] -= x[j]*a[i,j]
    x = aux

    print (x)
    gauss_seidel (a,x,b,n)


# gauss_jacobi(a,  [0,0,0,0,0], [1,1,1,5,2], 5)

print (x)
l, u = reduce_to_lu (x.astype(float),n)
print ("\n\n", l, "\n\n", u, "\n\n\n", np.matmul(l,u)==x )
print (np.matmul(l,u))

# r,r2 = [],[]
# for i in range(0,n):
#     r.append(solve_lower_tri (l,b[i],n))
# for i in range(0,n):
#     r2.append(solve_upper_tri (u,r[i],n))
# print (r2)
