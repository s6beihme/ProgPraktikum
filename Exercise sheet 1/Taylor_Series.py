# -*- coding: utf-8 -*-
"""
Created on Sat Oct 27 22:39:12 2018

@author: Benjamin Ihme
"""
from math import floor,factorial,exp,pi

# approximation of exponential function
# for given argument z go through the exponential series until N where 
# |z|<=1+N/2 and 2/(N+1)!*|z|^(N+1)

#NOTE: ONLY WORKS FOR ARGUMENTS NOT TOO BIG (THEN IT BECOMES SOMEWHAT INACCURATE)

def exp_approx(z):
    
    #check for type of z
    if type(z)!=int and type(z)!=float and type(z)!=complex:
        raise TypeError("Argument has to be int, float or complex")
    
    #set N such that |z|<=1+N/2
    N=2*floor(abs(z))
    
    #compute variables to possibly increase later
    N_plus_1_fact=factorial(N+1)
    z_to_N_plus_1=abs(z)**(N+1)
    
    # increase N until approximation of exponential will be accurate enough
    while 0.000001<(z_to_N_plus_1/N_plus_1_fact)*2:
        N_plus_1_fact*=(N+2)
        z_to_N_plus_1*=abs(z)
        N+=1
        
    #now approximate exp(z)
    i_fact=1
    z_to_i=1
    result=1
    
    for i in range(1,N+1):
        i_fact*=i
        z_to_i*=z
        result+=(1.0/i_fact)*z_to_i
    
    return result

#test 
print(exp(30)-exp_approx(30))
    
    