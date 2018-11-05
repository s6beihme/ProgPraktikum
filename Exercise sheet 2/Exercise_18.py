# -*- coding: utf-8 -*-
"""
Created on Mon Nov  5 21:45:16 2018

"""
import sympy as sp

# scalar product in C[-1,1] intended for polynomials
# first two arguments (f1 and f2) are sympy expressions, 
# the third argument (x) is the sympy Symbol with which both are defined 
# and with which they are to be integrated
def scalar_prod(f1,f2,x):
    #integrate using sympy integrate function (first argument is a sympy
    # 'function', second argument is tuple (symbol over which to integrate,
    # lower boundary, upper boundary))
    return sp.integrate(f1*f2, (x,-1,1))



# get input number n from user (afterwards compute first n legendre polynomials)
# (cast input to int)
n=int(input("how many legendre polynomials do you want?\n"))

x=sp.Symbol('x')

# L0=P0=1, (sum from 1 to (1-1)=0)=0 => L1=P1=x
L=[1]

for i in range(n-1):
    L.append(L[-1]*x)
    
    
# NOTE: GIVEN FORMULA IS NOT CORRECT, SUM HAS TO GO FROM 0 TO m-1 !!    
# now convert entries to legendre polynomials using corrected formula
for m in range(1,n):
    L[m]=(L[m]-sum([(scalar_prod(L[m],L[i],x)/scalar_prod(L[i],L[i],x))*L[i] for i in range(0,m)])).simplify()

print("the ")
print("the first", n, " legendre polynomials according to the corrected formula are:")
print(L)


#now compute legendre polynomials using three-term-recursion (faster for big n)
print("\nnow using the three-term-recursion (m+1)Lm=(2m+1)*x*L(m-1)-m*L(m-2):")
L2=[1,x]

for m in range(2,n):
    L2.append((sp.S(2*m-1)/sp.S(m)*x*L2[m-1]-(sp.S(m-1)/sp.S(m))*L2[m-2]).simplify())
print("\nthe first", n, "legendre polynomials are:")
print(L2)
print("note that three-term-recursion gives scaled results!")
print(scalar_prod(L2[0],L2[2],x))
