# -*- coding: utf-8 -*-
"""
Created on Sun Oct 21 20:24:36 2018

@author: Benjamin Ihme
"""

#TODO: LU DECOMPOSITIN AND INVERSE

#NOTE: I dont use __getitem__ and __setitem__ because that would increase runtime and typing "entries"
# a few times is not so bad
class Matrix:
    def __init__(self, entries=[[0]]):
        
        # check for mistakes and if necessary set the matrix to predefined value 
        if type(entries)!= list:
            if type(entries)==int or type(entries)==float or type(entries)==complex:
                self.entries=[[entries]]
            else:
                raise TypeError("incompatible type of entry for matrix")
        elif len(entries)==1 and len(entries[0])==0: #if given matrix is empty
            self.entries=[[0]] #set it to zero
        elif len(entries)==0: #if given matrix is empty
            self.entries=[[0]] #set it to zero
        elif type(entries[0])!=list: #check for 
            if type(entries[0])==int or type(entries[0])==float or type(entries[0])==complex:
                self.entries=[[x] for x in entries]
            else:
                raise TypeError("incompatible type of entry for matrix")
        else:
            k=len(entries[0]) #check if all rows have the same length
            for i in range(len(entries)):
                if len(entries[i])!= k:
                    raise Exception("all rows of matrix have to have same length")
                for j in range(len(entries[i])):
                    if type(entries[i][j])!=int and type(entries[i][j])!=float and type(entries[i][j])!=complex:
                        raise TypeError("only int, float and complex are accepted types  for entries of matrix")
            
            self.entries=entries
        
    def height(self):
        return len(self.entries)
    def width(self):
        return len(self.entries[0])
    
    #method to print matrix
    def __str__(self):
        s=''
        for i in range(self.height()):
            for j in range(self.width()):
                s= s+str(self.entries[i][j])+', '
            s=s+'\n'
        return s
        
    # method to multiply two matrices
    def __mul__(self,other):
        #check for compatible type
        if type(other)!= Matrix:
            raise TypeError("matrix can only be multiplied with other matrix")
        
        #check for compatible dimensions
        if self.width()!=other.height():
            raise Exception("Incompatible Dimensions for matrix multiplication")
        
        #multiply matrices
        new_entries=[[0 for x in range(other.width())] for n in range(self.height())] #0-entries of right dimension
        for i in range(self.height()):
            for k in range(other.width()):
                for j in range(self.width()):
                    new_entries[i][k]=new_entries[i][k]+self.entries[i][j]*other.entries[j][k]
        return Matrix(new_entries)
    
    # method to compare two matrices for equality
    def __eq__(self,other):
        
        #first check for equality despite other data type
        if type(other)!= Matrix:
            if type(other)==int or type(other)==float or type(other)==complex:
                if self.height()==1 and self.width()==1:
                    return self.entries[0][0]==other
                
            elif type(other)==list:
                if len(other)==0: return false
                elif len(other)==1:
                    return self.entries[0][0]==other[0]
                
                #check wether other is a "vector" equal to our matrix
                elif self.width()==1 and self.height()==len(other): 
                    for i in range(self.height()):
                        if self.entries[i][0]!=other[i]:
                            return False
                    return True
                
            else: return False
        else: return self.entries==other.entries
        
        # method to add row multiplier*r2 to row r1
    def add_row_to_row(self, r1, r2, multiplier):
        #check for false rows
        if r2<1 or r1<1 or r2>self.height() or r1>self.height():
            raise Exception("trying to add nonexistent row")
            
        #add mlutiplier*r2 to r1
        self.entries[r1-1]=[self.entries[r1-1][i]+multiplier*self.entries[r2-1][i] for i in range(self.width())]
    
    # method to multiply row r ba multiplier
    def mult_row(self,r, multiplier):
        #check for false row
        if r<1 or r>self.height():
            raise Exception("trying to multiply nonexistent row")
        
        #multiply row with multiplier
        self.entries[r-1]=[multiplier*self.entries[r-1][i] for i in range(self.width())]
        
    # method to swap row r1 with row r2
    def swap_rows(self,r1,r2):
        #check for false rows
        if r2<1 or r1<1 or r2>self.height() or r1>self.height():
            raise Exception("trying to swap nonexistent row")
        
        #swap rows r1 and r2
        if r1!=r2: # prevents useles computation
            self.entries[r1-1],self.entries[r2-1]=self.entries[r2-1],self.entries[r1-1]
    
    # method to compute inverse of square matrix
    def inverse(self):
        #check for correct dimensions
        if self.height()!=self.width():
            raise Exception("trying to find inverse of non square matrix")
        
        d=self.height() #to be used multiple times
        
        #build 1-Matrix of same dimension
        one_coef=[[1 if i==j else 0 for i in range(d)] for j in range(d)]
        E=Matrix(one_coef)
        
        #make copy of self
        copy=Matrix(self.entries)
        
        #apply steps used to form copy to the 1-Matrix to the 1 matrix
        for i in range(d):
            if copy.entries[i][i]!=0:
                
                E.mult_row(i+1, 1/copy.entries[i][i])
                copy.mult_row(i+1, 1/copy.entries[i][i])
              
                for j in list(range(i))+list(range(i+1, d)): #go through all rows except for i
                    E.add_row_to_row(j+1, i+1, (-1.0)*copy.entries[j][i]) #do same action on E
                    copy.add_row_to_row(j+1, i+1, (-1.0)*copy.entries[j][i]) #multiply so entry in i'th column is 0
                    
            else:
                nonzero_found=False
                nonzero_line_index=0
                for j in range(i+1, d): #go through all rows except for i
                    if copy.entries[j][i]!=0: #if non-zero entry is found below
                        nonzero_found=True
                        nonzero_line_index=j+1
                        break
                if nonzero_found==False: #if nonzero entry wasnt fount matrix isnt invertible
                    raise Exception("trying to invert non-invertable matrix")
                else: #if nonzero entry was found
                    
                    E.swap_rows(nonzero_line_index, i+1) #also swap in E
                    copy.swap_rows(nonzero_line_index, i+1) #swap the rows so nonzero is on diagonal
                    
                    E.mult_row(i+1, 1/copy.entries[i][i])
                    copy.mult_row(i+1, 1/copy.entries[i][i])
                    
                    for j in list(range(i))+list(range(i+1, d)): #go through all rows except for i
                        E.add_row_to_row(j+1, i+1, (-1.0)*copy.entries[j][i]) #do same action on E
                        copy.add_row_to_row(j+1, i+1, (-1.0)*copy.entries[j][i]) #multiply so entry in i'th column is 0
                        
        return E
                    
    #method to find LU Decomposition of square matrix
    def LU_Decomposition(self):
        
        #first check for if matrix is square
        if self.height()!=self.width():
            raise Exception("LU Decomposition only for square matrix")
            
        # create copy of matrix to perform actions on. it will be transformed 
        #to the Upper Matrix and store the entries of lower Matrix in Lower 
        #part, as that part of upper Matrix becomes 0
        U=Matrix(self.entries)
        
        d=self.height() #to be used multiple times
        
        #variables to help by searching
        row_max_entry=0
        max_entry=0
        
        # perform Decomposition
        for i in range(d):
            
            row_max_entry=i+1 #set var back to 0 so you can look for next biggest abs
            max_entry=abs(U.entries[i][i]) #set variable to current diagonal
            
            #go through all rows beneath row i+1 and search for maximal absolute entry
            for j in range(i+1,d):
                if abs(U.entries[j][i])>max_entry:
                    max_entry=abs(U.entries[j][i])
                    row_max_entry=j+1
            
            #check if max_entry is 0
            if max_entry==0:
                raise ValueError("trying to find LU Decomposition of non-regular matrix")
            
            #swap row with maximal absolute entry onto diagonal
            if row_max_entry!=i+1: 
                U.swap_rows(i+1, row_max_entry)
                
            #go through rows below and do special row to row addition, 
            #so that lower Matrix is not effected
            for j in range(i+1, d):
                
                #save entry of row j+1 in column i+1 because it will be changed,
                #end we need it later to add to Lower matrix
                b=U.entries[j][i]
                
                #go through all columns from i+1 on and subtract i+1'th row 
                #so that entry in i+1'th column becomes 0
                for k in range(i,d):
                    U.entries[j][k]=U.entries[j][k]+((-1.0)*b/U.entries[i][i])*U.entries[i][k]
                
                U.entries[j][i]=1.0*b/U.entries[i][i]
        
        return U

#Bsp LU Decomposition (be aware that the permutations of original matrix 
#are not represented in the output, so there is some permutation of rows P
#such that P*A=L*U) 
    
# =============================================================================
# M2=Matrix([[3,2,1],
#            [6,6,3],
#            [9,10,6]])  
# 
# print(M2.LU_Decomposition())
# 
# M=Matrix([[0,1],
#           [1,0]])
# print(M.LU_Decomposition())
# =============================================================================
                    
            
            
                
            
                        
            
            
#Example for inverse Matrix
#m=Matrix([[3 ,-1,2 ],
#          [-3,4 ,-1],
#          [-6,5,-2]])
#print(m.inverse())


        

#Example for assotiative matrix multiplication        
        
#m1=Matrix([[1,2,3,4]])

#m2=Matrix([[1,2  ],
#           [2,3  ],
#           [4,6  ],
#           [7,1.3]])
#    
#m3=Matrix([[3,5,7,3,1],
#           [6,2,9,0,0]])

#print((m1*m2)*m3)
#print(m1*(m2*m3))