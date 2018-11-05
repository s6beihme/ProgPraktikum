# -*- coding: utf-8 -*-
"""
Created on Sun Nov  4 16:11:24 2018

"""
import sympy as sp
e=sp.Symbol("e")

A=sp.Matrix([[1+e*sp.cos(2/e), -e*sp.sin(2/e)],
             [-e*sp.sin(2/e), 1+e*sp.cos(2/e)]])

#eigenvects is a function in sympy. it returns a list of tuples
#each tuple contains an eigenvalue, its algebraic multiplicity and the
#corresponding eigenvoctor(s)
l=A.eigenvects()
l1=l[0][0]
l2=l[1][0]
phi1=l[0][2]
phi2=l[1][2]
print("the eigenvalues are:\nl1=",l1,",\nl2=",l2)
#print and simplify eigenvectors
print("the eigenvectors are:\nphi1=",phi1, "\n=",sp.Matrix([[phi1[0][0].simplify()],[phi1[0][1].simplify()]]))
print("phi2=",phi2, "\n=", sp.Matrix([[phi2[0][0].simplify()],[phi2[0][1].simplify()]]))
# limit(f(x),x,t) returnes lim(x->t) f(x)
# use that function fo find A(e), li(e) and phii(e) for e->0
print("lim(e->0) e*cos(2/e)=", sp.limit(e*sp.cos(2/e),e,0))
print("=", )
print("lim(e->0) e*sin(2/e)=", sp.limit(e*sp.sin(2/e),e,0))
print("=>")
print("lim(e->0) A(e)=", sp.Matrix([[1+sp.limit(e*sp.cos(2/e),e,0), -sp.limit(e*sp.sin(2/e),e,0)],[-sp.limit(e*sp.sin(2/e),e,0), 1+sp.limit(e*sp.cos(2/e),e,0)]]))
print("lim(e->0) l1(e)=", sp.limit(l1,e,0))
print("lim(e->0) l2(e)=", sp.limit(l2,e,0))
print("phi1 and phi2 are constant")
