# -*- coding: utf-8 -*-
"""
Created on Sun Nov  4 18:43:28 2018

"""

import sympy as sp

u=sp.Function('u')
x=sp.Symbol('x')
# There is a bit of a problem: dsolve doesnt seem to care about boundary conditions
# (ics) and returnes the result with the constant C1, which is a 'module' object,
# as python tells you if you execute the following line
#C1=sp.symbol('C1')
# therefore you can not solve for C1 without some extra steps (see line 24)
print("the solution to u'(x)+u(x)=x, u(0)=1 is:")
#use dsolve from sympy to solve 
a=sp.dsolve(u(x).diff(x)+u(x)-x, ics={u(x).subs(x,0):1}).simplify()

#a is an object of type 'Equality' and the two sides of  the equality can be 
# called with 'args'
print(a.args[1])

#now find C1 manually so that u(0)=1
print("C1=",sp.solveset(a.args[1].subs(x,0)-1),"\n")

print("the solution to u''(x)=u(x), u(0)=0, u'(1)=1 is:")
a=sp.dsolve(sp.diff(u(x),x,x)-u(x), ics={u(x).subs(x,0):0, u(x).diff(x).subs(x,1):1})
print(a.args[1])
#compute C1 and C2 even more manually
c1,c2=sp.symbols('c1,c2')
print("if we set x=0 in the solution we get:")
print(a.args[1].subs(x,0))
print("if we set x=1 in the derivative of the solution we get:")
print(sp.diff(a.args[1],x).subs(x,1))

#now solve c1+c1=0 and -c1*exp(-1)+exp(1)*c2 as a linear system of equations
M=sp.Matrix([[1,1],
             [-sp.exp(-1), sp.exp(1)]])
b=sp.Matrix([0,-1])
c_sol=sp.linsolve((M,b),[c1,c2])
print("solving with the BVP we get:")
print("=> C1=",list(c_sol)[0][0], "\nC2=",list(c_sol)[0][1])
