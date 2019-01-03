# -*- coding: utf-8 -*-
"""
Created on Wed Jan  2 15:32:51 2019

@author: Benjamin Ihme
"""

#a)

import numpy as np
from math import sqrt
import os
import argparse

#get paths
parser=argparse.ArgumentParser(description="Driven-Cavity-Problem")
parser.add_argument('--input', type=str, required=True, help="Input file containing the parameters for the driven cavity problem")
parser.add_argument('--output', type=str, required=True, help="Base name for output files that result from driven cavity problem")
args=parser.parse_args()
inputname=args.input
outputname=args.output

#write all the functions from task a)


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
def apply_boundary_conditions_UV_2(imax, jmax, U, V):
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
        
#applies boundary conditions V and to U according to task b)
def apply_boundary_conditions_UV(imax, jmax, U, V):
    for j in range(1,jmax+1):
        U[0][j]=0
        U[imax][j]=0
        V[0][j]=-V[1][j]
        V[imax+1][j]=-V[imax][j]
    for i in range(1,imax+1):
        U[i][0]=-U[i][1]
        U[i][jmax+1]=2-U[i][jmax]
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
#has no return value
def calculate_F_G(imax, jmax, delt, delx, dely, alpha, GX, GY, U, V, F, G):
    #apply boundary conditions first
    for j in range(1, jmax+1):
        F[0][j]=U[0][j]
        F[imax][j]=U[imax][j]
    for i in range(1, imax):
        G[i][0]=V[i][0]
        G[i][jmax]=V[i][jmax]
    #now the rest
    for i in range(1,imax):
        for j in range(1, jmax+1):
            ddudxx=(U[i+1][j]-2*U[i][j]+U[i-1][j])/(delx*delx)
            ddudyy=(U[i][j+1]-2*U[i][j]+U[i][j-1])/(dely*dely)
            duudx=(1/delx)*((((U[i][j]+U[i+1][j])/2)**2)-(((U[i-1][j]+U[i][j])/2)**2))+alpha*(1/(4*delx))*((abs(U[i][j]+U[i+1][j])*(U[i][j]-U[i+1][j]))-(abs(U[i-1][j]+U[i][j])*(U[i-1][j]-U[i][j])))                            
            duvdy=(1/(dely*4))*((V[i][j]+V[i+1][j])*(U[i][j]+U[i][j+1])-(V[i][j-1]+V[i+1][j-1])*(U[i][j-1]+U[i][j]))+alpha*(1/(dely*4))*((abs(V[i][j]+V[i+1][j])*(U[i][j]-U[i][j+1]))-(abs(V[i][j-1]+V[i+1][j-1])*(U[i][j-1]-U[i][j])))
            F[i][j]=U[i][j]+delt*(((1/Re)*(ddudxx+ddudyy))-duudx-duvdy+GX)
    for i in range(1, imax+1):
        for j in range(1, jmax):
            duvdx=(1/(delx*4))*((U[i][j]+U[i][j+1])*(V[i][j]+V[i+1][j])-(U[i-1][j]+U[i-1][j+1])*(V[i-1][j]+V[i][j]))+alpha*(1/(delx*4))*((abs(U[i][j]+U[i][j+1])*(V[i][j]-V[i+1][j]))-(abs(U[i-1][j]+U[i-1][j+1])*(V[i-1][j]-V[i][j])))
            dvvdy=(1/dely)*((((V[i][j]+V[i][j+1])/2)**2)-(((V[i][j-1]+V[i][j])/2)**2))+alpha*(1/(4*dely))*((abs(V[i][j]+V[i][j+1])*(V[i][j]-V[i][j+1]))-(abs(V[i][j-1]+V[i][j])*(V[i][j-1]-V[i][j])))
            ddvdxx=(V[i+1][j]-2*V[i][j]+V[i-1][j])/(delx*delx)
            ddvdyy=(V[i][j+1]-2*V[i][j]+V[i][j-1])/(dely*dely)
            G[i][j]=V[i][j]+delt*(((1/Re)*(ddvdxx+ddvdyy))-duvdx-dvvdy+GY)


#cacluates the right hand side of formula (14)
#RHS matrix has to exist already
#has no return value
def calculate_RHS(imax, jmax, delt, delx, dely, F, G, RHS):
    for i in range(1, imax+1):
        for j in range(1, jmax+1):
            RHS[i][j]=(1/delt)*(((F[i][j]-F[i-1][j])/delx)+((G[i][j]-G[i][j-1])/dely))

#applies boundary conditios (19) to P
def apply_boundary_conditions_Pressure(imax, jmax, P):
    for j in range(1, jmax+1):
        P[0][j]=P[1][j]
        P[imax+1][j]=P[imax][j]
    for i in range(1, imax+1):
        P[i][0]=P[i][1]
        P[i][jmax+1]=P[i][jmax]
        
def make_B_to_a_copy_of_A(imax, jmax, B, A):
    for i in range(0, imax+2):
        for j in range(0, jmax+2):
            B[i][j]=A[i][j]
        
#calculates P using boundary conditions (19) and SOR (15)
#P has to exist already
#returns number of iterations and norm of residuum
def calculate_Pressure_with_SOR(imax, jmax, itermax, delx, dely, eps, omg, RHS, P):
    apply_boundary_conditions_Pressure(imax, jmax, P)
    #make second matrix to store the old values
    P_old=np.zeros([imax+2, jmax+2], float)
    make_B_to_a_copy_of_A(imax, jmax, P_old, P)
    
    it=0
    res_squared=eps*eps+1
    while it<itermax and res_squared>eps*eps:
        apply_boundary_conditions_Pressure(imax, jmax, P)
        apply_boundary_conditions_Pressure(imax, jmax, P_old)
        #do 1 SOR cycle
        for i in range(1, imax+1):
            for j in range(1, jmax+1):
                P[i][j]=(1-omg)*P_old[i][j]+(omg/(2*((1/delx**2)+(1/dely**2))))*(((P_old[i+1][j]+P[i-1][j])/delx**2)+((P_old[i][j+1]+P[i][j-1])/dely**2)-RHS[i][j])
        make_B_to_a_copy_of_A(imax, jmax, P_old, P)
        #calculate residuum (16)
        res_squared=sum(  [sum(  [(((P[i+1][j]-2*P[i][j]+P[i-1][j])/delx**2)+((P[i][j+1]-2*P[i][j]+P[i][j-1])/dely**2)-RHS[i][j])**2 for j in range(1, jmax+1)]  ) for i in range(1, imax+1)]  )/(imax*jmax)
        it+=1
# =============================================================================
#     if(res_squared>eps*eps):
#         print("SOR didnt yield a result that was close enough")
# =============================================================================
    return it, sqrt(res_squared)





#calculate velocities U and V according to (10) and (11)
def calculate_U_and_V(imax, jmax, delt, delx, dely, F, G, P, U, V):
    for i in range(1, imax):
        for j in range(1, jmax+1):
            U[i][j]=F[i][j]-(delt/delx)*(P[i+1][j]-P[i][j])
    for i in range(1, imax+1):
        for j in range(1, jmax):
            V[i][j]=G[i][j]-(delt/dely)*(P[i][j+1]-P[i][j])

def write_output_into_file(filename, xlength, ylength, imax, jmax, U, V, P, U2, V2):
    #colculate U2 and V2 as the values of U and V in the middle of cells
    for i in range(1, imax+1):
        for j in range(1, jmax+1):
            U2[i][j]=(U[i][j-1]+U[i][j])/2
            V2[i][j]=(V[i][j]+V[i+1][j])/2
    #write data into file
    with open(filename, 'w') as myfile:
        myfile.write(str(xlength))
        myfile.write("\n")
        myfile.write(str(ylength))
        myfile.write("\n")
        myfile.write(str(imax))
        myfile.write("\n")
        myfile.write(str(jmax))
        myfile.write("\n")
        for i in range(1,imax+1):
            for j in range(1, jmax+1):
                myfile.write(str(U2[i][j]))
                myfile.write(" ")
            myfile.write("\n")
        for i in range(1,imax+1):
            for j in range(1, jmax+1):
                myfile.write(str(V[i][j]))
                myfile.write(" ")
            myfile.write("\n")
        for i in range(1,imax+1):
            for j in range(1, jmax+1):
                myfile.write(str(P[i][j]))
                myfile.write(" ")
            myfile.write("\n")
        

#b)

#ALGORITHM:

#first create new directory to store resulting files in
#delete its content after using it, or ambiguities might occur


if not os.path.exists("Files"):
    os.makedirs("Files")

imax,jmax,xlength,ylength,delt,t_end,tau,del_vec,eps,omg,alpha,itermax,GX,GY,Re,UI,VI,PI=read_parameters_from_file(inputname)
del_vec2=del_vec

t=0
delx,dely=calculate_delx_dely(imax,jmax,xlength,ylength)
U,V,P=initialize_fields_UVP(imax,jmax, UI, VI, PI)
#matices that will be used in the algorithm
F=np.zeros([imax+2, jmax+2], float)
G=np.zeros([imax+2, jmax+2], float)
RHS=np.zeros([imax+2, jmax+2], float)
U2=np.zeros([imax+2, jmax+2], float)
V2=np.zeros([imax+2, jmax+2], float)
number_of_files=0
while t<t_end:
    delt=calculate_delt(tau,Re,delx,dely,U,V, delt)
    
    apply_boundary_conditions_UV(imax,jmax,U,V)
    #print(U)
    
    calculate_F_G(imax, jmax, delt, delx, dely, alpha, GX, GY, U, V, F, G)
    calculate_RHS(imax, jmax, delt, delx, dely, F, G, RHS)
    
    it, res=calculate_Pressure_with_SOR(imax, jmax, itermax, delx, dely, eps, omg, RHS, P)
    
    calculate_U_and_V(imax, jmax, delt, delx, dely, F, G, P, U, V)

    if(t>del_vec):
        number_of_files+=1
        write_output_into_file("Files/"+outputname+"_{0:03}".format(number_of_files), xlength, ylength, imax, jmax, U, V, P, U2, V2)
        del_vec=t+del_vec2
    t+=delt
number_of_files+=1
write_output_into_file("Files/"+outputname+"_{0:03}".format(number_of_files), xlength, ylength, imax, jmax, U, V, P, U2, V2)

