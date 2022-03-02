# -*- coding: utf-8 -*-
"""
Spyder Editor

This is a temporary script file.
"""
import numpy as np
import matplotlib.pyplot as plt
import math

class Laplas():
    def __init__(self,x, y, h, func_fi_0, func_fi_l ,func_psi_0 , func_psi_l,eps):
        self.eps=eps
        self.x= x
        self.x_0=0
        self.y=y
        self.y_0=0
        self.h=h
        self.func_fi_0=func_fi_0
        self.func_fi_l=func_fi_l
        self.func_psi_0=func_psi_0
        self.func_psi_l=func_psi_l
        self.X=np.arange(self.x_0,self.x+self.h, self.h)
        self.Y=np.arange(self.y_0,self.y+self.h,  self.h)
        self.len_X=len(self.X)
        self.len_Y=len(self.Y)
        print("init")
        
    def StartFillMatrix(self):
        T=np.zeros((self.len_X, self.len_Y))
        T[:,0]=self.func_psi_0(self.Y[::-1])
        T[:,-1]=self.func_psi_l(self.Y[::-1])
        T[-1,:]=self.func_fi_l(self.X[::-1])
        T[0,:]=self.func_fi_0(self.X[::-1])

        return T
    
    def iterative_procedures(self,w):
        T=self.StartFillMatrix()
        T_new=np.zeros((self.len_X, self.len_Y))
        T2=self.StartFillMatrix()
        k=0
        # for i in range(1,self.len_X-1):
#                 for j in range(1,self.len_Y-1):
#                     T[i,j]=(1/4)*(T[i-1,j]+T[i+1,j]+T[i,j-1]+T[i,j+1])
        

        while (True) :
            k=k+1
            delta=[]
            for i in range(1,self.len_X-1):
                for j in range(1,self.len_Y-1):
                    T[i,j]=(1/4)*(T[i-1,j]+T[i+1,j]+T[i,j-1]+T[i,j+1])
            for i in range(1,self.len_X-1):
                for j in range(1,self.len_Y-1):
                    T2[i,j]=(1/4)*(T[i-1,j]+T[i+1,j]+T[i,j-1]+T[i,j+1])

            for i in range(0,self.len_X):
                for j in range(0,self.len_Y):
                    T_new[i,j]=T[i,j]+w*(T2[i,j]-T[i,j])      
                    
            for i in range(0,self.len_X):
                for j in range(0,self.len_Y):
                    delta.append(abs(T_new[i,j]-T[i,j]))
            
            if(max(delta)<self.eps):
                print("break")
                break
            else:
                for i in range(0,self.len_X):
                    for j in range(0,self.len_Y):
                        T[i,j]=T_new[i,j]
               
        return [T_new, k]
    
x,h,y,eps=1,0.2,1,0.01

func_psi_0 = lambda y: 20*y
func_psi_l = lambda y: [30*math.cos((math.pi*i)/2) for i in y]
func_fi_0  = lambda x: 20*(x**2)
func_fi_l  = lambda x: [30*math.cos((math.pi*i)/2) for i in x]

Test=Laplas(x, y, h, func_fi_0, func_fi_l ,func_psi_0 , func_psi_l,eps)
mat=Test.StartFillMatrix()
mass_k=[]
mass_w=[]
for i in range(1,20):
#i=1
    m=i/10
    mass_w.append(m)
    print(m)
    
    _return=Test.iterative_procedures(m)
    res1=_return[0]
    k=_return[1]
    mass_k.append(k)
    print(k)
    if (i==19):
        X=np.linspace(0, Test.x, num = Test.len_X)
        Y=np.linspace(0, Test.y, num = Test.len_Y)
        X,Y = np.meshgrid(X,Y)
        fig, ax = plt.subplots()
        ax.set_title("Разрез по y")
        ax.plot(X, res1[:,1],label="y=0.2")
        ax.plot(X, res1[:,3],label="y=0.6")
        ax.plot(X, res1[:,5],label="y=1")
        ax.set_xlabel("y")
        ax.set_ylabel("U")
        ax.legend()
        fig, ax = plt.subplots()
        ax.set_title("Разрез по x")
        ax.plot(Y, res1[1,:],label="x=0.2")
        ax.plot(Y, res1[3,:],label="x=0.6")
        ax.plot(Y, res1[5,:],label="x=1")
        ax.set_xlabel("y")
        ax.set_ylabel("U")
    
        ax.legend()
        fig, ax = plt.subplots(subplot_kw={"projection": "3d"})
        Z=res1
        surf = ax.plot_surface(X, Y, res1, cmap='coolwarm',linewidth=0)
        ax.set_xlabel("x")
        ax.set_ylabel("y")
        ax.set_zlabel("U")
        ax.view_init(40, 10)
        plt.show()

fig, ax = plt.subplots()
ax.plot(mass_w, mass_k,label="")
ax.set_xlabel("w")
ax.set_ylabel("k")
plt.show()