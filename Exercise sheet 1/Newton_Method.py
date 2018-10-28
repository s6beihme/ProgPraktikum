# -*- coding: utf-8 -*-
"""
Created on Sat Oct 27 23:25:14 2018

@author: Benjamin Ihme
"""

# simple implementation of the newton method for one dimensional continuously
# differetiable functions, where the derivative has to be given explicitly
# to stop in case of not convergent sequence do 100 iterations at most

from math import sqrt

def newton(f, f_prime, x):
    #check for types
    if callable(f)==0 or callable(f_prime)==0:
        raise TypeError("first two arguments of newton fct have to be callable")
        
    if type(x)!=int and type(x)!=float:
        raise TypeError("third argument has to be int or float")
    
    #do first iteration
    next_x=x-(f(x)/f_prime(x))
    i=1
    
    #iterate while |x_(n+1)-x_n|>10^-6
    while i<=100 and abs(next_x-x)>0.000001:
        x=next_x
        next_x=x-(f(x)/f_prime(x))
        i+=1
        
    return next_x

def f(x):
    return x*x-2

def f_prime(x):
    return 2*x

print(newton(f,f_prime,1)-sqrt(2))
