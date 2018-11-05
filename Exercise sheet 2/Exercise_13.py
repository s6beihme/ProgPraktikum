# -*- coding: utf-8 -*-
"""
Created on Sun Nov  4 14:43:21 2018

"""
import sympy as sp



a=sp.Symbol("a")
A=sp.Matrix([[1,1,1],
             [1,0,1],
             [a,-1,1]])

betta=sp.Symbol("betta")
b=sp.Matrix([3,betta,3])


print("a)")
print("For Ax=0 to only have the trivial solution x=0, A has to be injective, and thus invertible.\n")
print("A is invertible <=> det(A)!=0\n")
print("det(A)=",A.det())
print("=> det(A)!=0 <=> a isnt in",sp.solve(a-1,a))

print("\nb)")

print("The column vectors of a square Matrix are linearly dependent if det(A)=0\n")
print("det(A) for a=0 is", sp.det(A.subs(a,0)))
print("=> the column vectors of A for a=0 are not linearly dependent for a=0")

print("\nc)")
print("if a is so that A is invertible, then x=A^(-1)*b is a solution for Ax=b for all betta.")
print("=> let a=1 (see a)).")
x,y,z=sp.symbols("x,y,z")
print("the equation Ax=b (for random betta) hat the set of solutions:",sp.linsolve((A.subs(a,1),b),[x,y,z]))
print("therefore Ax=b has no solution for a=1 and betta in the real numbers")

print("\nd)")
print("for a=-3 and betta=0 Ax=b has the set of solutions: ")
print(sp.linsolve((A.subs(a,-3),b.subs(betta,0)),[x,y,z]))

