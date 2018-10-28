# -*- coding: utf-8 -*-
"""
Created on Sun Oct 28 12:57:26 2018

@author: Benjamin Ihme
"""

from scipy import integrate
from math import sin,pi
import numpy as np

def f2(x):
    return 3**(3*x-1)
def f3(x):
    return exp(x**2)

# integral of sin with 3 differnt methods
print("quad:", integrate.quad(sin,0,pi))
print("romberg:", integrate.romberg(sin,0,pi))
x=np.linspace(0,pi,100)
print("trapezodial", integrate.trapz([sin(z) for z in x],x))

# integral of f2 with 3 differnt methods
print("\nquad:", integrate.quad(f2,0,2))
print("romberg:", integrate.romberg(f2,0,2))
x=np.linspace(0,2,100)
print("trapezodial", integrate.trapz([f2(z) for z in x],x))

# integral of f3 with 3 differnt methods
print("\nquad:", integrate.quad(f3,0,1))
print("romberg:", integrate.romberg(f3,0,1))
x=np.linspace(0,1,100)
print("trapezodial", integrate.trapz([f3(z) for z in x],x))