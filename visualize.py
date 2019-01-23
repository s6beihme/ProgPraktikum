import argparse
import numpy as np
import matplotlib.pyplot as plt


parser=argparse.ArgumentParser(description="Driven-Cavity-Problem")
parser.add_argument('--input', type=str, required=True, help="Input file containing the parameters for the driven cavity problem")
args=parser.parse_args()
inputname=args.input


def read_data_from_file(inputname):
    with open(inputname, 'r') as myfile:
        xlength=float(myfile.readline())
        ylength=float(myfile.readline())
        imax=int(myfile.readline())
        jmax=int(myfile.readline())
        U=[]
        for i in range(imax):
            line=myfile.readline().split(" ")
            del line[-1]
            U.append([float(x) for x in line])
        V=[]
        for i in range(imax):
            line=myfile.readline().split(" ")
            del line[-1]
            V.append([float(x) for x in line])
        P=[]
        for i in range(imax):
            line=myfile.readline().split(" ")
            del line[-1]
            P.append([float(x) for x in line])
    return xlength, ylength, imax, jmax, U, V, P

xlength, ylength, imax, jmax, U, V, P = read_data_from_file(inputname)
plt.quiver(V,U, headwidth=2 )
plt.matshow(P)
plt.show()