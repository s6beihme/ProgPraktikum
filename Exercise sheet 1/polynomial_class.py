# -*- coding: utf-8 -*-
"""
Created on Sun Oct 21 16:33:38 2018

@author: Benjamin Ihme
"""

# a simple polynomial class
class Polynomial:
    
    #Constructor where the coefficients are given in a list
    def __init__(self, coef2=[0]):
        if type(coef2)!= list:
            raise ValueError("coefficients have to be given in a list")
        elif len(coef2)==0:
            self.coef=[0]
        else:    
            self.coef=list(coef2)
        
        #now go through coefficients and if the highest one is 0, then delete it
        while len(self.coef)>1 and self.coef[-1]==0:
            self.coef.pop()
    
    # str method to print polynomials
    def __str__(self):
        #create a string representing the polynomial
        s=''
        for i in range(len(self.coef)-1, -1, -1): 
            s='+('+str(self.coef[i])+')x^'+str(i)+s
        return s
    
    # method to evaluate polynomial at given position
    def __call__(self,x):
        #check if it type of argument is correct
        if type(x)!= int and type(x)!=float and type(x)!= complex:
            raise TypeError("Can only compute polynomial at a complex point")
        
        r=0
        #add ai*x^i for all 0<=i<=grad(p)
        for i in range(len(self.coef)):
            r=r+self.coef[i]*(x**i)
        return r
    
    # method to add two polynomials
    def __add__(self, other):
        
        #first check for compatability
        if type(other) != Polynomial:
            #if type isnt Polynomial check for int, float and complex
            if type(other) == int or type(other)==float or type(other)==complex: 
                p=Polynomial(self.coef)
                p.coef[0]=p.coef[0]+other
                return p
            # if type is not compatible raise ValueError
            else:
                raise TypeError("Only Polynomials, ints, floats and comlex numbers can be added to a polynomial")
                
        # create new polynomial, with the added coefficients of self and other
        if len(self.coef)>=len(other.coef):
            new_coef=self.coef
            for i in range(len(other.coef)):
                new_coef[i]=new_coef[i]+other.coef[i]
        else:
            new_coef=other.coef
            for i in range(len(self.coef)):
                new_coef[i]=new_coef[i]+self.coef
        
        return Polynomial(new_coef)
    
    # method to multiply two polynomials
    def __mul__(self,other):
        #first check for compatability
        if type(other) != Polynomial:
            #if type isnt Polynomial check for int, float and complex
            if type(other) == int or type(other)==float or type(other)==complex: 
                new_coef=[other*a for a in self.coef]
                return Polynomial(new_coef)
            # if type is not compatible raise ValueError
            else:
                raise TypeError("Only Polynomials, ints, floats and comlex numbers can be multiplied with a polynomial")
                
        # create new polynomial
        new_coef=[0 for x in range(len(self.coef)+len(other.coef)-1)]
        for i in range(len(self.coef)): #go through all coeff. in self
            for j in range(len(other.coef)): #and in other
                new_coef[i+j]=new_coef[i+j]+self.coef[i]*other.coef[j] #ad the multiple of the two to the right position
        
        return Polynomial(new_coef)
        
    # method to check if two polynomials are equal
    def __eq__(self,other):
        #first check for compatability
        if type(other) != Polynomial:
            #if type isnt Polynomial check for int, float and complex
            if type(other) == int or type(other)==float or type(other)==complex: 
                if len(self.coef)==1 and self.coef[0]==other:
                    return True
            # if type is not compatible return False
            else:
                return False
        #now compare coeff.
        return self.coef==other.coef
        
    # method to find derivative of polynomial
    def derivative(self):
        new_coef=[]
        for i in range(len(self.coef)-1):
            new_coef.append((i+1)*self.coef[i+1])
        return Polynomial(new_coef)
    
    # method to find the antiderivative of a polynomial
    #assuming that the constant term is allway 0
    def antiderivative(self):
        new_coef=[0]
        for i in range(len(self.coef)):
            new_coef.append(self.coef[i]/(i+1.0))
        return Polynomial(new_coef)
        
    
#test derivative and antiderivative for p(x)=3x^2+2x+1
p1=Polynomial([1,2,3])
print(p1)
print(p1.derivative())
print(p1.antiderivative())

    
    
    
    
    
    
    
    
    
    

