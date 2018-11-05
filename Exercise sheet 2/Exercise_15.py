# -*- coding: utf-8 -*-
"""
Created on Sun Nov  4 17:17:53 2018

"""

import sympy as sp
import numpy as np
import matplotlib.pyplot as plt

x,a=sp.symbols("x,a")
expr=x**3+3*x-a

print("the set of possible solutions for x**3+3*x-a=0 for a being a real number are:\n")
#use the sympy function solve to solve expr=0 for variable x
print(sp.solveset(expr,x,domain=sp.S.Complexes))

print("\nfor any real a the equation only has 1 real solution, because")
print("d(x**3+3*x-a)/dx =", sp.diff(expr,x), ">0 for all real x\n")

print("A visual representation of the real solution for a given a:")


a_axis=np.linspace(-500,500,50)
x_axis=np.array([list(sp.solveset(expr.subs(a,i),x,domain=sp.S.Reals))[0] for i in a_axis])
plt.plot(a_axis,x_axis, label='real solution with respect to a')
plt.xlabel('a')
plt.ylabel('x')
plt.legend()
plt.show()

