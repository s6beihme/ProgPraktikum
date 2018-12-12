# -*- coding: utf-8 -*-
"""
Created on Sun Nov 25 19:21:30 2018

"""
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import axes3d
import numpy as np
def plot_from_file(filename):
    with open(filename, 'r') as f:
        size=int(f.readline())
        if size==2:
            x_axis=[]
            y_axis=[]
            for line in f:
                temp=line.partition(" ")
                x_axis.append(float(temp[0]))
                y_axis.append(float(temp[2]))
            plt.plot(x_axis, y_axis, label='B-Spline Curve')
            plt.legend()
            plt.show()
        if size==3:
            fig=plt.figure()
            ax=fig.add_subplot(111, projection='3d')
            
            X,Y,Z=[],[],[]
            for line in f:
                temp=line.partition(" ")
                temp2=temp[2].partition(" ")
                X.append(float(temp[0]))
                Y.append(float(temp2[0]))
                Z.append(float(temp2[2]))
            X=np.array(X)
            Y=np.array(Y)
            Z=np.array([Z])
            ax.plot_wireframe(X,Y,Z)
            plt.show()
    
plot_from_file("splinefile1.txt")