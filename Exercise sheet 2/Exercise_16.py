# -*- coding: utf-8 -*-
"""
Created on Sun Nov  4 18:12:34 2018

"""
import sympy as sp

# Newto method. first argument is an expression of sympy variables, 
# second argument is the variable of f, the third is a real number 
def newton(f, x, x0):
    #check for right type of x
    if type(x0)!=int and type(x0)!=float:
        raise TypeError("third argument has to be int or float")
    
    #use diff from sympy to compute the derivative of f
    f_prime=sp.diff(f,x)
    
    #do first iteration
    # convert result to float because otherwise expression might "explode"
    next_x0=float(x0-(f.subs(x,x0)/f_prime.subs(x,x0)))
    i=1
    #iterate while |x_(n+1)-x_n|>10^-6
    while i<=100 and abs(next_x0-x0)>0.000001:
        x0=next_x0
        next_x0=float(x0-(f.subs(x,x0)/f_prime.subs(x,x0)))
        i+=1
        
    return next_x0

x=sp.Symbol('x')
f=sp.exp(x)+2*x
a=newton(f,x,1)
print("the solution to exp(x)+2x is:")
print(float(a))

f2=sp.cosh(x)-2*x

a=newton(f2,x,2)
print("\none solution to cosh(x)-2x is:")
print(a)
a=newton(f2,x,-2)
print("\nthe other one is:")
print(a)


#TODO: b)
