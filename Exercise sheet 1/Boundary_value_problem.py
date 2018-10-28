# -*- coding: utf-8 -*-
"""
Created on Sun Oct 28 16:04:46 2018

@author: Benjamin Ihme
"""
import numpy as np
from scipy.integrate import solve_bvp

# Solve the boundary value problem u''(x)=-|u(x)|

# model of u''
def model(x,u):
    #unpack y
    u0=u[0]
    u1=u[1]
    
    du0dx=u1
    du1dx=-abs(u0)
    return np.array([du0dx, du1dx])

# boundary conditions: bc(ua,ub+2) = 0
def bc(ua,ub):
    return np.array([ua[0], ub[0]+2])

#initial mesh
x=np.linspace(0,4,2)

# here store initial guesses u(0)=u[0,0], u'(0)=u[1,0]
u=np.zeros((2,x.size)) 
u[1,0]=1 #u'(0)>0

#solve bvp
res=solve_bvp(model, bc, x,u)

#plot solution
x_plot=np.linspace(0,4,50)
y_plot=res.sol(x_plot)[0]
plt.plot(x_plot, y_plot, 'r-', label='solution for u\'(0)>0')
plt.legend()
plt.show()

u=np.zeros((2,x.size)) 
u[1,0]=-1 #u'(0)<0

#solve bvp
res=solve_bvp(model, bc, x,u)

#plot solution
x_plot=np.linspace(0,4,50)
y_plot=res.sol(x_plot)[0]
plt.plot(x_plot, y_plot, 'b-', label='solution for u\'(0)<0')
plt.legend()
plt.show()
