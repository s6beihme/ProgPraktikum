# -*- coding: utf-8 -*-
"""
Created on Wed Jan  2 15:32:51 2019

"""
#reads parameters from file (file has to be in the format from exercise sheet)
#returns imax,jmax,xlength,ylength,delt,t_end,tau,del_vec,eps,omg,alpha,itermax,GX,GY,Re,UI,VI,PI
#which are all integers of floats 
#UI,VI,PI are constants with which the matrices are filled
def read_parameters_from_file(filename):
    with open(filename, 'r') as myfile:
        imax=int(str(myfile.readline()).partition("=")[2])
        jmax=int(str(myfile.readline()).partition("=")[2])
        xlength=float(str(myfile.readline()).partition("=")[2])
        ylength=float(str(myfile.readline()).partition("=")[2])
        delt=float(str(myfile.readline()).partition("=")[2])
        t_end=float(str(myfile.readline()).partition("=")[2])
        tau=float(str(myfile.readline()).partition("=")[2])
        del_vec=float(str(myfile.readline()).partition("=")[2])
        eps=float(str(myfile.readline()).partition("=")[2])
        omg=float(str(myfile.readline()).partition("=")[2])
        alpha=float(str(myfile.readline()).partition("=")[2])
        itermax=int(str(myfile.readline()).partition("=")[2])
        GX=float(str(myfile.readline()).partition("=")[2])
        GY=float(str(myfile.readline()).partition("=")[2])
        Re=int(str(myfile.readline()).partition("=")[2])
        UI=float(str(myfile.readline()).partition("=")[2])
        VI=float(str(myfile.readline()).partition("=")[2])
        PI=float(str(myfile.readline()).partition("=")[2])
        return imax,jmax,xlength,ylength,delt,t_end,tau,del_vec,eps,omg,alpha,itermax,GX,GY,Re,UI,VI,PI
print(type(read_parameters_from_file("examplefile.txt")[2]))
    