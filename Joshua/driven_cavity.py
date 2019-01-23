#!/usr/bin/env python
# coding: utf-8

# In[1]:


import numpy as np
import os
import argparse

#parsen
parser=argparse.ArgumentParser(description="Driven-Cavity-Problem")
parser.add_argument('--input', type=str, required=True, help="Input file containing the parameters for the driven cavity problem")
parser.add_argument('--output', type=str, required=True, help="Base name for output files that result from driven cavity problem")
args=parser.parse_args()
inputname=args.input
outputname=args.output

#max value in matrix
def Umax(U,imax,jmax):
    Umax=0.0
    for i in range (0,imax+2):
        for j in range (0,jmax+2):
            if(abs(U[i][j])>Umax):
                Umax=abs(U[i][j])
    return Umax
            
def Vmax(V,imax,jmax):
    Vmax=0.0
    for i in range (0,imax+2):
        for j in range (0,jmax+2):
            if(abs(V[i][j])>Vmax):
                Vmax=abs(V[i][j])
    return Vmax


#derivatives 
def _u2_x(U,i,j,delx,alpha):
    return 1/delx*(((U[i][j]+U[i+1][j])/2)**2-((U[i-1][j]+U[i][j])/2)**2)+alpha*1/delx*(((abs(U[i][j]+U[i+1][j]))/2)*((U[i][j]-U[i+1][j])/2)-((abs(U[i-1][j]+U[i][j]))/2)*((U[i-1][j]-U[i][j])/2))
                                                                                        
def _uv_y(U,V,i,j,dely,alpha):
    return 1/dely*((V[i][j]+V[i+1][j])/2*(U[i][j]+U[i][j+1])/2-(V[i][j-1]+V[i+1][j-1])/2*(U[i][j-1]+U[i][j])/2)+alpha*1/dely*((abs(V[i][j]+V[i+1][j]))/2*(U[i][j]-U[i][j+1])/2-(abs(V[i][j-1]+V[i+1][j-1]))/2*(U[i][j-1]-U[i][j])/2) 
                                                                                    
def _2u_x2(U,i,j,delx):
    return (U[i+1][j]-2*U[i][j]+U[i-1][j])/(delx**2)
    
def _2u_y2(U,i,j,dely):
    return (U[i][j+1]-2*U[i][j]+U[i][j-1])/(dely**2)
    
def _p_x(P,i,j,delx):
    return (P[i+1][j]-P[i][j])/delx
    
def _v2_y(V,i,j,dely,alpha):
    return 1/dely*(((V[i][j]+V[i][j+1])/2)**2-((V[i][j-1]+V[i][j])/2)**2)+alpha*1/dely*(((abs(V[i][j]+V[i][j+1]))/2)*((V[i][j]-V[i][j+1])/2)-((abs(V[i][j-1]+V[i][j]))/2)*((V[i][j-1]-V[i][j])/2))
    
def _uv_x(U,V,i,j,delx,alpha):
    return 1/dely*((U[i][j]+U[i][j+1])/2*(V[i][j]+V[i+1][j])/2-(U[i-1][j]+U[i-1][j+1])/2*(V[i-1][j]+V[i][j])/2)+alpha*1/dely*((abs(U[i][j]+U[i][j+1]))/2*(V[i][j]-V[i+1][j])/2-(abs(U[i-1][j]+U[i-1][j+1]))/2*(V[i-1][j]-V[i][j])/2)
                                                                                       
def _2v_x2(V,i,j,delx):
    return (V[i+1][j]-2*V[i][j]+V[i-1][j])/(delx**2)
    
def _2v_y2(V,i,j,dely):
    return (V[i][j+1]-2*V[i][j]+V[i][j-1])/(dely**2)
    
def _p_y(P,i,j,dely):
    return (P[i][j+1]-P[i][j])/dely
     
def output(filename, xlength, ylength, imax, jmax, U, V, P):
    #compute centered version of staggered grid
    Ucentered=np.zeros((imax+2,jmax+2))
    Vcentered=np.zeros((imax+2,jmax+2))
    for i in range(1, imax+1):
        for j in range(1, jmax+1):
            Ucentered[i][j]=(U[i][j-1]+U[i][j])/2
            Vcentered[i][j]=(V[i][j]+V[i+1][j])/2
    #create output file
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
                myfile.write(str(Ucentered[i][j]))
                myfile.write("\n")
        for i in range(1,imax+1):
            for j in range(1, jmax+1):
                myfile.write(str(Vcentered[i][j]))
                myfile.write("\n")
        for i in range(1,imax+1):
            for j in range(1, jmax+1):
                myfile.write(str(P[i][j]))
                myfile.write("\n")

if not os.path.exists("Files"):
    os.makedirs("Files")
                
                
def _input(filename):
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

imax,jmax,xlength,ylength,delt,t_end,tau,del_vec,eps,omg,alpha,itermax,GX,GY,Re,UI,VI,PI=_input(inputname)

#init values and fields
t=0.0
tvis=del_vec
delx=xlength/imax
dely=ylength/jmax
filenumber=0

U=np.zeros((imax+2,jmax+2))
V=np.zeros((imax+2,jmax+2))
P=np.zeros((imax+2,jmax+2))
for i in range (0,imax+2):
    for j in range (0,jmax+2):
        U[i][j]=UI
        V[i][j]=VI
        P[i][j]=PI
RHS=np.zeros((imax+2,jmax+2))
F=np.zeros((imax+2,jmax+2))
G=np.zeros((imax+2,jmax+2))



print("start")
        
while t < t_end:
    #compute delt
    if(tau>=0):
        u_max=Umax(U,imax,jmax)
        if(u_max==0):
            u_max=delx/(Re/(2*(1/delx**2+1/dely**2)))
        v_max=Vmax(V,imax,jmax)
        if(v_max==0):
            v_max=delx/(Re/(2*(1/delx**2+1/dely**2)))
        delt=min(Re/(2*(1/delx**2+1/dely**2)),delx/u_max,dely/v_max)   
    
    #Boundary Values U,V
    for j in range (1,jmax+1):
        V[0][j]=-V[1][j]
        V[imax+1][j]=-V[imax][j]
    for i in range (1,imax+1):
        U[i][0]=-U[i][1]
        U[i][jmax+1]=-U[i][jmax]
    for i in range (1,imax+1):
        U[i][jmax+1]=2.0-U[i][jmax]
        
    #copmute F, G
    for i in range (0,imax):
        for j in range (0,jmax+1):
            F[i][j]=U[i][j]+delt*(1/Re*(_2u_x2(U,i,j,delx)+_2u_y2(U,i,j,dely))-_u2_x(U,i,j,delx,alpha)-_uv_y(U,V,i,j,dely,alpha)+GX)
           
    for i in range (0,imax+1):
        for j in range (0,jmax):
             G[i][j]=V[i][j]+delt*(1/Re*(_2v_x2(V,i,j,delx)+_2v_y2(V,i,j,dely))-_uv_x(U,V,i,j,delx,alpha)-_v2_y(V,i,j,dely,alpha)+GY)
    
    #Boundary Values F,G
            
    for j in range (1,jmax+1):
        F[imax][j]=U[imax][j]
        F[0][j]=U[0][j]
    
    for i in range (1,imax+1):
        G[i][jmax]=V[i][jmax]
        G[i][0]=V[i][0]
     
    #compute RHS
    for i in range (1,imax+1):
        for j in range (1,jmax+1):
            RHS[i][j]=1/delt*((F[i][j]-F[i-1][j])/delx+(G[i][j]-G[i][j-1])/dely)
    
    #compute P
    it=0
    res=1.0
    while ((it<itermax) and (res>eps)):
        
        #Pressure boundary values
        P[0][j]=P[1][j]
        for j in range (1,jmax+1):
            P[imax+1][j]=P[imax][j]
        P[i][0]=P[i][1]
        for i in range (1,imax+1):
            P[i][jmax+1]=P[i][jmax]

        #SOR-Cycle
        P_old=np.zeros((imax+2,jmax+2))
        for i in range (0,imax+2):
            for j in range (0,jmax+2):
                P_old[i][j]=P[i][j]
        for i in range (1,imax+1):
            for j in range (1,jmax+1):
                P[i][j]=(1-omg)*P_old[i][j]+omg/(2*(1/(delx**2)+1/(dely**2)))*((P_old[i+1][j]+P[i-1][j])/(delx**2)+(P_old[i][j+1]+P[i][j-1])/(dely**2)-RHS[i][j])
        for i in range (1,imax+1):
            for j in range (1,jmax+1):
                res=((P[i+1][j]-2*P[i][j]+P[i][j-1])/(delx**2)+(P[i][j+1]-2*P[i][j]+P[i][j-1])/(dely**2)-RHS[i][j])**2/(imax*jmax)
        res=res**(1/2)
        it=it+1
    
    #compute U,V    
    for i in range (1,imax):
        for j in range (1,jmax+1):
            U[i][j]=F[i][j]-delt/delx*(P[i+1][j]-P[i][j])
    for i in range (1,imax+1):
        for j in range (1,jmax):
            V[i][j]=G[i][j]-delt/dely*(P[i][j+1]-P[i][j])           
    
    #create file
    print(t)
    t=t+delt
    if(t>=tvis):
        filenumber=filenumber+1
        output("Files/"+outputname+"_{0:03}".format(filenumber), xlength, ylength, imax, jmax, U, V, P)
        tvis=t+del_vec
    
    
print("ende")  


# In[ ]:





# In[ ]:




