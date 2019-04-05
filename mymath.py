"""
muller method
"""

import cmath
def muller(f,x0,x1,x2,tol=1e-09,nmax=50):
    h1 = x1 - x0
    h2 = x2 - x1
    t1 = (f(x1)-f(x0))/h1
    t2 = (f(x2)-f(x1))/h2
    d = (t2-t1)/(h2+h1)
    #for i in range(nmax):
    i = 0
    while i<=nmax:
        b = t2 + h2 * d
        t = cmath.sqrt(b**2 - 4 * f(x2) * d)
        if abs(b-t)<abs(b+t):
            h = -2.0*f(x2)/(b+t)
        else:
            h = -2.0*f(x2)/(b-t)
        p = x2+h
        if abs(h)<tol:
            return i,p
        x0 = x1
        x1 = x2
        x2 = p
        h1 = x1 - x0
        h2 = x2 - x1
        t1 = (f(x1)-f(x0))/h1
        t2 = (f(x2)-f(x1))/h2
        d = (t2-t1)/(h2+h1)
        i+=1
    print("Over maximum times")
    return i,p

def func(x):
    return 16.0*x**4 - 40.0*x**3+5.0*x**2+20.0*x+6.0

"""
x = muller(func,0.5+0j,-0.5+0j,0+0j)
print(x)
"""

