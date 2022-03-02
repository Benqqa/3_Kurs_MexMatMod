# -*- coding: utf-8 -*-
"""
Spyder Editor

This is a temporary script file.
"""
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import cm
from mpl_toolkits.mplot3d import Axes3D
import math
import pandas as pd

class Struna():
    def __init__(self, x, t, h, dt, func_x_0, func2_x_0 , d_func_x_0 , func_0_t , func_1_t):
        self.x= x
        self.x_0=0
        self.t=t
        self.t_0=0
        self.h=h
        self.func_x_0=func_x_0
        self.func2_x_0=func2_x_0
        self.d_func_x_0=d_func_x_0
        self.func_0_t=func_0_t
        self.func_1_t=func_1_t
        self.X=np.arange(self.x_0,self.x, self.h)
        self.dt=self.t/(len(self.X))
        self.Time=np.arange(self.t_0,self.t,  self.dt)
        self.len_X=len(self.X)
        self.len_T=len(self.Time)
        print("init")
    def StartFillMatrix(self):
        T=np.zeros((self.len_X, self.len_T))
        T[:,0]=self.func_x_0(self.X)
        T[0,:]=self.func_0_t(self.Time)
        T[-1,:]=self.func_1_t(self.Time)
        
        #T[:,1]
        for i in range(1,self.len_X-1):
            print(i)
            print(self.X[i])
            T[i,1]=T[i,0]+self.d_func_x_0(self.X[i])*self.dt+((self.dt**2)/(2*(self.h**2)))*(T[i,0]-2*T[i,0]+T[i,0])

        return T
    def explicit_schema(self):
        T=self.StartFillMatrix()
        c=1
        for k in range(1,self.len_T-1):
            for i in range(1,self.len_X-1):
                T[i,k+1]=(((c**2)*(self.dt**2))/(self.h**2))*(T[i+1,k]-2*T[i,k]+T[i-1,k])+2*T[i,k]-T[i,k-1]
                
        return T
    def implicit_schema(self):
        T=self.StartFillMatrix()
        A=1/(self.h**2)
        C=1/(self.h**2)
        B=(2*(self.dt**2)+(self.h**2))/((self.h**2)*(self.dt**2))
        F=np.zeros((self.len_X))
        P=np.zeros((self.len_X))
        Q=np.zeros((self.len_X))
        for k in range(1,self.len_T-1):
            for i in range(0, self.len_X):
                F[i]=((2*T[i,k])/(self.dt**2))-(T[i,k-1]/(self.dt**2))
            P[0]=C/B
            Q[0]=F[0]/B
            for i in range(1,self.len_X):
                P[i]=C/(B-A*P[i-1])
                Q[i]=(F[i]+A*Q[i-1])/(B-A*P[i-1])
            for i in range(self.len_X-2,0,-1):
                T[i,k+1]=P[i]*T[i+1,k+1]+Q[i]    
                
        return T
x,h,t,dt=1,0.1,0.5,0.01
func_x_0  = lambda x: [math.cos((math.pi*i)/2) for i in x]
func2_x_0 = lambda x: math.cos((math.pi*x)/2)
d_func_x_0= lambda x: x**2
func_0_t  = lambda t: 1+2*t
func_1_t  = lambda t: 0
Test=Struna(x, t, h, dt, func_x_0, d_func_x_0 ,func2_x_0 , func_0_t , func_1_t)
res1=Test.explicit_schema()
res2=Test.implicit_schema()
print ("explicit_schema \n",pd.DataFrame(res1))
print ("implicit_schema \n",pd.DataFrame(res2))
X=np.linspace(0, Test.x, num = Test.len_X)
Y=np.linspace(0, Test.t, num = Test.len_T)
#X,Y = np.meshgrid(X,Y)
X1=np.arange(0,x, h)
dt=t/(len(X1))
Y1=np.arange(0,t,  dt)
print('X',X)
print('Y',Y)
print('X1',X1)
print('Y1',Y1)
fig, ax = plt.subplots()
ax.plot(X, res1[:,2],label="t2="+str(Y1[2]))
ax.plot(X, res1[:,5],label="t5="+str(Y1[5]))
ax.plot(X, res1[:,8],label="t8="+str(Y1[8]))
ax.set_xlabel("x")
ax.set_ylabel("U")
ax.set_title("Разрезы в разные моменты времени(явная схема)")
ax.legend()
fig, ax = plt.subplots()
ax.plot(X, res2[:,2],label="t2="+str(Y1[2]))
ax.plot(X, res2[:,5],label="t5="+str(Y1[5]))
ax.plot(X, res2[:,8],label="t8="+str(Y1[8]))
ax.set_xlabel("x")
ax.set_ylabel("U")
ax.set_title("Разрезы в разные моменты времени(неявная схема)")
ax.legend()

fig, ax = plt.subplots()
ax.plot(X, res1[:,5],label="неявная схема")
ax.plot(X, res2[:,5],label="явная схема")
ax.set_xlabel("x")
ax.set_ylabel("U")
ax.set_title("Разрезы в один момент времени для явной и неявной схемы")
ax.legend()
fig, ax = plt.subplots()
ax.plot(X, res1[:,2],label="неявная схема")
ax.plot(X, res2[:,2],label="явная схема")
ax.set_xlabel("x")
ax.set_ylabel("U")
ax.set_title("Разрезы в один момент времени для явной и неявной схемы")
ax.legend()
fig, ax = plt.subplots(subplot_kw={"projection": "3d"})
Z=res1
surf = ax.plot_surface(X, Y, res1, cmap='coolwarm',linewidth=0)
ax.set_xlabel("x")
ax.set_ylabel("t")
ax.set_zlabel("U")
ax.set_title("явная схема")
ax.view_init(20, 40)
#ax.view_init(40, 20)

fig, ax = plt.subplots(subplot_kw={"projection": "3d"})
Z=res2
surf = ax.plot_surface(X, Y, res2, cmap='coolwarm',linewidth=0)
ax.set_xlabel("x")
ax.set_ylabel("t")
ax.set_zlabel("U")
ax.set_title("неявная схема")
ax.view_init(20, 40)
#ax.view_init(40, 20)
plt.show()









