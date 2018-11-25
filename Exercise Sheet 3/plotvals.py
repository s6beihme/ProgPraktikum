# -*- coding: utf-8 -*-
"""
Created on Sun Nov 25 19:21:30 2018

"""
import matplotlib.pyplot as plt
def plot_from_file(filename):
    x_axis=[]
    y_axis=[]
    with open(filename, 'r') as f:
        size=int(f.readline())
        for line in f:
            temp=line.partition(" ")
            x_axis.append(float(temp[0]))
            y_axis.append(float(temp[2]))
    plt.plot(x_axis, y_axis, label='B-Spline Curve')
    plt.legend()
    plt.show()
    
plot_from_file("C:/Users/ASUS/Documents/Studium/Prog Praktikum/Exercise Sheet 3/testfile2.txt")