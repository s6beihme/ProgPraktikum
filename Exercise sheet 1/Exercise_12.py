# -*- coding: utf-8 -*-
"""
Created on Sun Oct 28 17:26:26 2018

@author: Benjamin Ihme
"""
#


from numpy import *

n=51
x=linspace(0,1,n)
y=linspace(0,1,n)
L=zeros((n,n))

L[0,:]= y*(y-1)



iter=70
for _ in range(iter):
    for i in range(1,n-1):
        for j in range(1,n-1):
            Lt=L
            L[i,j]=(Lt[i+1,j]+Lt[i-1,j]+Lt[i,j+1]+Lt[i,j-1])/4
    #print(L)

#note: initially distribution is not symmetric to middle column of matrix
#However, with more iterations it becomes more and more symmetric, because
# the "advantage" of the entries to the right of the middle column becomes smaller
# with the "input" adding up in each iteration
print(L[1,1], L[1,3])
plt.matshow(L)

# I dont really know what you mean by plotting the matrix over the
# square determined by x and y. scaling down the axes by n or something?

#The program computes how something flows through some object that has static 
#boundaries and a static amount of input from the top.
# Each position in the object is effected only by its neighbors and in each
#iteration receives the mean value of them.
# at some point (after enough iterations) values reach a state of equilibrium

