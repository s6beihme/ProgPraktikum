# -*- coding: utf-8 -*-
"""
Created on Wed Jan  2 15:32:51 2019

"""
import numpy as np

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


#initializes U,V and P according to UI,VI,PI 
def initialize_fields_UVP(imax,jmax, UI, VI, PI):
    U=np.zeros([imax+2,jmax+2],float)
    V=np.zeros([imax+2,jmax+2],float)
    P=np.zeros([imax+2,jmax+2],float)
    for i in range(1,imax+1):
        for j in range(1,jmax+1):
            U[i][j]=UI
            V[i][j]=VI
            P[i][j]=PI
    return U,V,P


#applies boundary conditions (17, 18) to U and V
def apply_boundary_conditions_UV(imax, jmax, U, V):
    for j in range(1,jmax+1):
        U[0][j]=0
        U[imax][j]=0
        V[0][j]=-V[1][j]
        V[imax+1][j]=-V[imax][j]
    for i in range(1,imax+1):
        U[i][0]=-U[i][1]
        U[i][jmax+1]=-U[i][jmax]
        V[i][0]=0
        V[i][jmax]=0

def calculate_delx_dely(imax,jmax,xlength,ylength):
    return xlength/imax, ylength/jmax

def calculate_delt(tau,Re,delx,dely,U,V, delt):
    if(tau<0):
        return delt
    umax=float(max(np.amax(U), abs(np.amin(U))))
    vmax=float(max(np.amax(U), abs(np.amin(U))))
    if umax==0:
        if vmax==0:
            return tau*(Re/2.0)*(1.0/((1.0/(delx*delx))+(1.0/(dely*dely))))
        else:
            return tau*min((Re/2.0)*(1.0/((1.0/(delx*delx))+(1.0/(dely*dely)))), dely/vmax)
    else:
        if vmax==0:
            return tau*min((Re/2.0)*(1.0/((1.0/(delx*delx))+(1.0/(dely*dely)))), delx/umax)
        else:
            return tau*min((Re/2.0)*(1.0/((1.0/(delx*delx))+(1.0/(dely*dely)))), delx/umax, dely/vmax)

#calculates F and G from U and V and applies boundary conditions
#F and G have to exist already
def calculate_F_G(imax, jmax, delt, delx, dely, alpha, GX, GY, U, V, F, G):
    #apply boundary conditions first
    for j in range(1, jmax+1):
        F[0][j]=U[0][j]
        F[imax][j]=U[imax][j]
    for i in range(1, imax):
        G[i][0]=V[i][0]
        G[i][jmax]=V[i][jmax]
    #now the rest
    for i in range(1,imax+1):
        for j in range(1, jmax):
            ddudxx=(U[i+1][j]-2*U[i][j]+U[i-1][j])/(delx*delx)
            ddudyy=(U[i][j+1]-2*U[i][j]+U[i][j-1])/(dely*dely)
            duudx=(1/delx)*((((U[i][j]+U[i+1][j])/2)**2)-(((U[i-1][j]+U[i][j])/2)**2))+alpha*(1/(4*delx))*((abs(U[i][j]+U[i+1][j])*(U[i][j]-U[i+1][j]))-(abs(U[i-1][j]+U[i][j])*(U[i-1][j]-U[i][j])))                            
            duvdy=(1/(dely*4))*((V[i][j]+V[i+1][j])*(U[i][j]+U[i][j+1])-(V[i][j-1]+V[i+1][j-1])*(U[i][j-1]+U[i][j]))+alpha*(1/(dely*4))*((abs(V[i][j]+V[i+1][j])*(U[i][j]-U[i][j+1]))-(abs(V[i][j-1]+V[i+1][j-1])*(U[i][j-1]-U[i][j])))
            
            duvdx=(1/(delx*4))*((U[i][j]+U[i][j+1])*(V[i][j]+V[i+1][j])-(U[i-1][j]+U[i-1][j+1])*(V[i-1][j]+V[i][j]))+alpha*(1/(delx*4))*((abs(U[i][j]+U[i][j+1])*(V[i][j]-V[i+1][j]))-(abs(U[i-1][j]+U[i-1][j+1])*(V[i-1][j]-V[i][j])))
            dvvdy=(1/dely)*((((V[i][j]+V[i][j+1])/2)**2)-(((V[i][j-1]+V[i][j])/2)**2))+alpha*(1/(4*dely))*((abs(V[i][j]+V[i][j+1])*(V[i][j]-V[i][j+1]))-(abs(V[i][j-1]+V[i][j])*(V[i][j-1]-V[i][j])))
            ddvdxx=(V[i+1][j]-2*V[i][j]+V[i-1][j])/(delx*delx)
            ddvdyy=(V[i][j+1]-2*V[i][j]+V[i][j-1])/(dely*dely)
            
            F[i][j]=U[i][j]+delt*(((1/Re)*(ddudxx+ddudyy))-duudx-duvdy+GX)
            G[i][j]=V[i][j]+delt*(((1/Re)*(ddvdxx+ddvdyy))-duvdx-dvvdy+GY)
    return F,G


imax,jmax,xlength,ylength,delt,t_end,tau,del_vec,eps,omg,alpha,itermax,GX,GY,Re,UI,VI,PI=read_parameters_from_file("examplefile.txt")
U,V,P=initialize_fields_UVP(imax,jmax, UI, VI, PI)
apply_boundary_conditions_UV(imax,jmax,U,V)
delx,dely=calculate_delx_dely(imax,jmax,xlength,ylength)
delt=calculate_delt(tau,Re,delx,dely,U,V, delt)
F=np.zeros([imax+2, jmax+2], float)
G=np.zeros([imax+2, jmax+2], float)
F,G=calculate_F_G(imax, jmax, delt, delx, dely, alpha, GX, GY, U, V, F, G)
print(F)