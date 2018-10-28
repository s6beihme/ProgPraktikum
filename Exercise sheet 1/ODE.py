# -*- coding: utf-8 -*-
"""
Created on Sun Oct 28 13:43:51 2018

@author: Benjamin Ihme
"""

import numpy as np
import matplotlib.pyplot  as plt
from scipy.integrate import odeint


# function that returns u'(x) for given u(x) and x
print("first ODE:\n")

def model1(u,x):
    return x-u


x=np.linspace(0,10,40)
y=odeint(model1,1,x)
y2=[z-1+2*exp(-z) for z in x]  

#plot numeric result
plt.plot(x,y,'r-',label='approximation')
plt.legend()
plt.show()

#plot analytic result
plt.plot(x,y2,'b-',label='exact')
plt.legend()
plt.show()


#now solve second ODE by converting it so system of ODEs

print("\nsecond ODE:\n")


def model2(z,x):
    #unpack y
    u0=z[0]
    u1=z[1]
    
    du0dx=u1
    du1dx=((3*x+2)*u1+(6*x-8)*u0)/(3*x-1)
    return [du0dx, du1dx]

x2=np.linspace(0,10,40)
z0=[2,3]

#solve ODE using ODE-int
z=odeint(model2, z0,x2)
u=z[:,0]

#plot numeric result
plt.plot(x,u,'r-',label='approximation')
plt.legend()
plt.show()

#plot analytic result
u=[2*exp(2*t)-t*exp(-t) for t in x]
plt.plot(x,u,'b-', label='exact')
plt.legend()
plt.show()





    
    