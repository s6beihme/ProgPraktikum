# -*- coding: utf-8 -*-
"""
Created on Thu Oct 25 22:26:53 2018

@author: Benjamin Ihme
"""
# NOT A GENERAL FUNCTION TO USE NEWTONS METHOD! JUST ONE EXPLICIT COMPUTATION

import scipy
import scipy.linalg

# function of which we want to find the zero
def f(xyz):
    #checking if argument is list of length 3
    if type(xyz)!=scipy.ndarray:
        raise TypeError("Arument has to be list")
        
    if len(xyz)!=3:
        raise ValueError("List has to have length 3")
    
    #checking if types in list are correct (long)
    if type(xyz[0])!=scipy.int32 and type(xyz[0])!=scipy.float64 and type(xyz[0])!=scipy.complex128:
        raise TypeError("arguments in list have to be int, float or complex")
    if type(xyz[1])!=scipy.int32 and type(xyz[1])!=scipy.float64 and type(xyz[1])!=scipy.complex128:
        raise TypeError("arguments in list have to be int, float or complex")
    if type(xyz[0])!=scipy.int32 and type(xyz[1])!=scipy.float64 and type(xyz[2])!=scipy.complex128:
        raise TypeError("arguments in list have to be int, float or complex")
    
    # return the function value
    return scipy.array([9*xyz[0]**2+36*xyz[1]**2+4*xyz[2]**2-36,
                        xyz[0]**2-2*xyz[1]**2-20*xyz[2],
                        xyz[0]**2-xyz[1]**2+xyz[2]**2])
    
    
    
#the jacobian matrix of the function f
def J(xyz):
    #checking if argument is list of length 3
    if type(xyz)!=scipy.ndarray:
        raise TypeError("Arument has to be list")
        
    if len(xyz)!=3:
        raise ValueError("List has to have length 3")
    
    #checking if types in list are correct (long)
    if type(xyz[0])!=scipy.int32 and type(xyz[0])!=scipy.float64 and type(xyz[0])!=scipy.complex128:
        raise TypeError("arguments in list have to be int, float or complex")
    if type(xyz[1])!=scipy.int32 and type(xyz[1])!=scipy.float64 and type(xyz[1])!=scipy.complex128:
        raise TypeError("arguments in list have to be int, float or complex")
    if type(xyz[0])!=scipy.int32 and type(xyz[1])!=scipy.float64 and type(xyz[2])!=scipy.complex128:
        raise TypeError("arguments in list have to be int, float or complex")
        
    #return the jacobian matrix
    return scipy.array([[18*xyz[0],72*xyz[1],8*xyz[2]],
                        [2*xyz[0] ,-4*xyz[1], -20],
                        [2*xyz[0] ,-2*xyz[1], 2*xyz[2]]])
#print(f(scipy.array([1,1,0])))
    
#now approximate solution to equation f(v)=0
u0=scipy.array([1,1,0])
for i in range(10):
    u0=u0-scipy.linalg.inv(J(u0))@f(u0).T
    
print(u0,f(u0))

u0=scipy.array([-1,-1,0])
for i in range(10):
    u0=u0-scipy.linalg.inv(J(u0))@f(u0).T
    
print(u0,f(u0))

        
    