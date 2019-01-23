#!/usr/bin/env python
# coding: utf-8

# In[ ]:


import argparse
import numpy as np
import matplotlib.pyplot as plt

#parsen
parser=argparse.ArgumentParser(description="Driven-Cavity-Problem")
parser.add_argument('--input', type=str, required=True, help="Input file containing the parameters for the driven cavity problem")
args=parser.parse_args()
inputname=args.input

#reads data from file that has been created by driven_cavity.py
#returns xlenght, ylength, imax, jmax, U, V, P
def read_data_from_file(inputname):
    with open(inputname, 'r') as myfile:
        xlength=float(myfile.readline())
        ylength=float(myfile.readline())
        imax=int(myfile.readline())
        jmax=int(myfile.readline())
        U=np.zeros((imax,jmax))
        V=np.zeros((imax,jmax))
        P=np.zeros((imax,jmax))
        for i in range(imax):
            for j in range (imax):
                U[i][j]=myfile.readline()
        for i in range(imax):
            for j in range (imax):
                V[i][j]=myfile.readline()
        for i in range(imax):
            for j in range (imax):
                P[i][j]=myfile.readline()
    return xlength, ylength, imax, jmax, U, V, P

xlength, ylength, imax, jmax, U, V, P=read_data_from_file(inputname)
plt.quiver(V,U, headwidth=2 )
plt.matshow(P)
plt.show()
        



