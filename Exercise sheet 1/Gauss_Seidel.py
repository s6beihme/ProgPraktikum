# -*- coding: utf-8 -*-
"""
Created on Wed Oct 24 14:47:06 2018

@author: ASUS
"""

import scipy
import scipy.linalg

#Gauss-seidel method to approximate solutin for LSE Ax=b (a and b given as matrix)
def gauss_seidel(A,b,x0):
    #check for false datatype
    if type(A)!=scipy.ndarray or type(b)!=scipy.ndarray or type(x0)!=scipy.ndarray:
        raise TypeError("Arguments of gauss_seidel have to be given as scipy arrays")
    
    #check for false dimensions
    if len(A)!=len(A[0]) or len(A)!=len(b) or len(x0)!=len(A):
        raise Exception("Height of Matrix and vector have to be the same")
        
    #Create L, U and D
    D=scipy.diag([A[i][i] for i in range(len(A))])
    U=scipy.array([[A[i][j] if j>i else 0 for j in range(len(A))] for i in range(len(A))])
    L=scipy.array([[A[i][j] if j<i else 0 for j in range(len(A))] for i in range(len(A))])
    
    #compute (D+L)^-1 to be used multiple times later
    I=scipy.linalg.inv(D+L)
    
    # do recursive approximation
    for i in range(2):
        x0=I@(b-U@x0)
    
    return x0

#Example with given matrix
A=scipy.array([[4,3,0],
               [3,4,-1],
               [0,-1,4]])
b=scipy.array([[24,30,-24]]).T
x0=scipy.array([[3,3,3]]).T
print(len(A),len(b),len(x0))

print(gauss_seidel(A,b,x0))
    
    
    
    