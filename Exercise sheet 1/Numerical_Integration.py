# -*- coding: utf-8 -*-
"""
Created on Thu Oct 25 20:00:27 2018

@author: Benjamin Ihme
"""
from math import sin,pi,exp

# function to approximate the integral over a function f from a to be unsing
# the composite trapezodial rule with n subintervalls
def trapeze(f, a, b, n):
    #check wether f is callable / f is a function
    if not callable(f):
        raise TypeError("function given to 'trapeze' isnt callable")
    
    #check if a<b
    if a>b:
        raise ValueError("lower boundary has to be lower than upper boundary")
    
    #check if n is an integer and >=1
    if type(n)!=int:
        raise TypeError("last argument has to be an integer")
    if n<1:
        raise ValueError("last argument has to be >=1")
        
    #now do the computation
    h=(b-a)/n
    return h*(1/2*(f(a)+f(b))+sum([f(a+i*h) for i in range(1,n)]))

# define function that approximates the integral a given function from a to b 
# up to 6 decimal places, given the function, the boundaries and the exact 
# integral
def approx(f,a,b,ex_res):
    r=trapeze(f,a,b,1)
    n=1
    
    #reapproximate integral with one more interval as long as the approximation
    #isnt exact to 6 decimal places
    while abs(r-ex_res)>0.000001: 
        n+=1
        r=trapeze(f,a,b,n)
    return r

#example
#print(approx(sin,0,pi,2))

# improvement of function trapeze by giving it the result of previous function call
# cutting each interval in half and thus only compute function values at new
# points. Because this is only to be used by doubling number of intervalls each
# time it is called, I do not check, wether n is 2^m
def trapeze2 (f,a,b,n,prev_val):
    #check wether f is callable / f is a function
    if not callable(f):
        raise TypeError("function given to 'trapeze' isnt callable")
    
    #check if a<b
    if a>b:
        raise ValueError("lower boundary has to be lower than upper boundary")
    
    #check if n is an integer and >=1
    if type(n)!=int:
        raise TypeError("last argument has to be an integer")
    if n<1:
        raise ValueError("last argument has to be >=1")
        
    #now do the computation
    #to get the result that trapeze would return divide the previous value
    #by two (becouse interval size has been divided by two) and add 
    #h*(sum of values of f at all the new points)
    h=(b-a)/n
    return prev_val/2.0+h*sum([f(a+i*h) for i in range(1,n,2)])

#improvement of function approx by using trapeze2 (and trapeze for first step)
def approx2(f,a,b,ex_res):
    r=trapeze(f,a,b,1)
    n=1
    
    #reapproximate integral with double the number of intervals
    #as long as the approximation isnt exact to 6 decimal places
    while abs(r-ex_res)>0.000001: 
        n=n*2
        r=trapeze2(f,a,b,n,r)
    return r

#example
#print(approx2(sin,0,pi,2))

#store approximations for n=2^m for 1 <= m <= 10
# =============================================================================
# value_store=[trapeze(sin,0,pi,2)]
# n=2
# 
# for i in range(9):
#     n=n*2
#     value_store.append(trapeze2(sin,0,pi,n,value_store[-1]))
#     
# #print the results
# print("values of computations:")
# for i in range(len(value_store)):
#     print("\n",value_store[i])
# print("\ndeviations from exact result:")
# for i in range(len(value_store)):
#     print("\n",abs(value_store[i]-2))
# =============================================================================
    
# approximating with 2^n intervalls leads to an exponential convergence 
# towards the exact integral for n increasing (each step the deviation from
#the exact result is divided by roughly 4)
# approximatin with n intervalls leads to a quadratic convergence 
#towards the exact integral for n increasing

#check assumption for 3^(3x-1)
# =============================================================================
# value_store2=[trapeze(lambda x:3**(3*x-1),0,2,2)]
# n=2
# 
# for i in range(9):
#     n=n*2
#     value_store2.append(trapeze2(lambda x:3**(3*x-1),0,2,n,value_store2[-1]))
#     
# #print the results
# # exact value approximated with wolfram alpha
# print("values of computations:")
# for i in range(len(value_store2)):
#     print("\n",value_store2[i])
# print("\ndeviations from exact result:")
# for i in range(len(value_store2)):
#     print("\n",abs(value_store2[i]-73.628239664926))
# =============================================================================
    

# example: integral of e^(x^2) from 0 to 1
# by applying formula for derivation: 1155 intervalls are naccesary to be exact 
# to 6 decimal places
# =============================================================================
# print(trapeze(lambda x:exp(x*x), 0,1,1155))
# print(abs(trapeze(lambda x:exp(x*x), 0,1,1155)-1.46265174590718160))
# =============================================================================
    
