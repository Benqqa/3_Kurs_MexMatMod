# -*- coding: utf-8 -*-
"""
Spyder Editor

This is a temporary script file.
"""
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import cm
from matplotlib import colors
from mpl_toolkits.mplot3d import Axes3D
import math
import pandas as pd
from matplotlib.animation import FuncAnimation

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
        self.dt=dt
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
    def explicit_schema(self,ro,E,mu,S,m,M,nu,k1,k2):
        T=self.StartFillMatrix()
        c=1
        for k in range(1,self.len_T-1):
            for i in range(1,self.len_X-1):
                A=((self.dt**2)*E/(i*(self.h**2)))
                # T[i,k+1]=(((c**2)*(self.dt**2))/(self.h**2))*(T[i+1,k]-2*T[i,k]+T[i-1,k])+2*T[i,k]-T[i,k-1]
                T[i,k+1]=1/(ro+self.dt*E*mu) * ((A*i)*T[i-1,k]+(A*i)*T[i+1,k]-T[i,k]*((A*i)-2*ro+self.dt*E*mu)-ro*T[i,k-1])
                # B=E*S/(self.h**2)
                # C=M+m
                # D=k1+k2*nu
                # T[i,k+1]=1/(C/(dt**2)+D/dt)*(B*T[i+1,k]+B*T[i-1,k]+(D/self.dt+2*C/(self.dt**2)-2*B)*T[i,k]-C/(self.dt**2)*T[i,k-1])
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

x,h,t,dt=1,0.1,8,0.025
# func_x_0  = lambda x: [math.cos((math.pi*i)/2) for i in x]
# func2_x_0 = lambda x: math.cos((math.pi*x)/2)
# d_func_x_0= lambda x: x**2
# func_0_t  = lambda t: 0
# func_1_t  = lambda t: 0

func_x_0  = lambda x: 0.1
func2_x_0 = lambda x: 0.1
d_func_x_0= lambda x: 0
func_0_t  = lambda t: 0
func_1_t  = lambda t: 0
Test=Struna(x, t, h, dt, func_x_0, d_func_x_0 ,func2_x_0 , func_0_t , func_1_t)
res1=Test.explicit_schema(7000,51,5,0.002827,0.01,1,0.00015,0.015,1720)
# res2=Test.implicit_schema()
print ("explicit_schema \n",pd.DataFrame(res1))
# print ("implicit_schema \n",pd.DataFrame(res2))
X=np.linspace(0, Test.x, num = Test.len_X)
Y=np.linspace(0, Test.t, num = Test.len_T)
# X,Y = np.meshgrid(X,Y)
X1=np.arange(0,x, h)
dt=t/(len(X1))
Y1=np.arange(0,t,  dt)
# print('X',X)
# print('Y',Y)
# print('X1',X1)
# print('Y1',Y1)
fig, ax = plt.subplots()
for i in range(len(Y1)):
    ax.plot(X, res1[:,i],label="t"+str(i)+"="+str(Y1[i]))
    #,label="t"+str(i)+"="+str(Y1[i])
    

# ax.plot(X, res1[:,5],label="t5="+str(Y1[5]))
# ax.plot(X, res1[:,8],label="t8="+str(Y1[8]))
# ax.set(ylim=(-1000000, 1000000))
ax.set_xlabel("x")
ax.set_ylabel("U")
ax.set_title("?????????????? ?? ???????????? ?????????????? ??????????????(?????????? ??????????)")
ax.legend()
# fig, ax = plt.subplots()
# ax.plot(X, res2[:,2],label="t2="+str(Y1[2]))
# ax.plot(X, res2[:,5],label="t5="+str(Y1[5]))
# ax.plot(X, res2[:,8],label="t8="+str(Y1[8]))
# ax.set_xlabel("x")
# ax.set_ylabel("U")
# ax.set_title("?????????????? ?? ???????????? ?????????????? ??????????????(?????????????? ??????????)")
# ax.legend()
############################3
# fig, ax = plt.subplots()
# ax.plot(X, res1[:,5],label="?????????? ??????????")
# # ax.plot(X, res2[:,5],label="?????????? ??????????")
# ax.set_xlabel("x")
# ax.set_ylabel("U")
# ax.set_title("?????????????? ?? ???????? ???????????? ?????????????? ?????? ?????????? ?? ?????????????? ??????????")
# ax.legend()
##########################
# fig, ax = plt.subplots()
# ax.plot(X, res1[:,2],label="?????????? ??????????")
# # ax.plot(X, res2[:,2],label="?????????? ??????????")
# ax.set_xlabel("x")
# ax.set_ylabel("U")
# ax.set_title("?????????????? ?? ???????? ???????????? ?????????????? ?????? ?????????? ?? ?????????????? ??????????")
# ax.legend()
# fig, ax = plt.subplots(subplot_kw={"projection": "3d"})
# Z=res1
# surf = ax.plot_surface(X, Y, res1, cmap='coolwarm',linewidth=0)
# ax.set_xlabel("x")
# ax.set_ylabel("t")
# ax.set_zlabel("U")
# ax.set_title("?????????? ??????????")
# ax.view_init(20, 40)
#ax.view_init(40, 20)

# fig, ax = plt.subplots(subplot_kw={"projection": "3d"})
# Z=res2
# surf = ax.plot_surface(X, Y, res2, cmap='coolwarm',linewidth=0)
# ax.set_xlabel("x")
# ax.set_ylabel("t")
# ax.set_zlabel("U")
# ax.set_title("?????????????? ??????????")
# ax.view_init(20, 40)
#ax.view_init(40, 20)
plt.show()
###########################
print(len(Y))
Frames=len(Y)
fig = plt.figure()
ax = plt.axes(xlim=(0, 4), ylim=(-2, 2))
line, = ax.plot([], [], lw=3)
 
def init():
    line.set_data([], [])
    return line,
def animate(i):
    line.set_data(X, res1[:,i])
    return line,
 
anim = FuncAnimation(fig, animate, init_func=init,
                                frames=Frames, interval=20, blit=True )
 
 
anim.save('radius.gif', writer='imagemagick')
###########################3
# fig = plt.figure()
# ax = fig.add_subplot(projection='3d')
# # Create the mesh in polar coordinates and compute corresponding Z.
# r=np.linspace(0, Test.x, num = Test.len_X)
# p = np.linspace(0,2*np.pi, num = Test.len_X)
# R, P = np.meshgrid(r, p)
# Z = res1[:,1]*R

# # Express the mesh in the cartesian system.
# X, Y = R*np.cos(P), R*np.sin(P)

# # Plot the surface.
# ax.plot_surface(X, Y, Z, cmap=plt.cm.winter)

# # Tweak the limits and add latex math labels.
# ax.set_xlim(0, 4)
# ax.set_zlim(0, 4)
# ax.set_ylim(0, 2)

# ax.set_xlabel("x")
# ax.set_ylabel("y")
# ax.set_zlabel("z")
# plt.show()

# fig = plt.figure()
# ax = fig.add_subplot(projection='3d')
# surf = ax.plot_surface([], [], np.array([[],[]]),cmap=plt.cm.winter)
# ax.set_xlim(0, 4)
# ax.set_zlim(0, 4)
# ax.set_ylim(0, 2)
# def initZ():
#     surf.set_data([], [], np.array([[],[]]))
#     return surf
# def animateZ(i):
#     # Create the mesh in polar coordinates and compute corresponding Z.
#     r=np.linspace(0, Test.x, num = Test.len_X)
#     p = np.linspace(0,2*np.pi, num = Test.len_X)
#     R, P = np.meshgrid(r, p)
#     Z = res1[:,i]*R
#     print(Z)
#     # Express the mesh in the cartesian system.
#     X, Y = R*np.cos(P), R*np.sin(P)
    
#     surf.set_data(X, Y, Z)
#     return surf
 
# anim = FuncAnimation(fig, animateZ, init_func=initZ,
#                                 frames=len(Y), interval=200, blit=True)
 
 
# anim.save('sine_waaaaave.gif', writer='imagemagick')
###################################
cmap = plt.cm.brg
# plt.cm.winter
# colors.ListedColormap(['red', 'blue'])
bounds = np.linspace(-0.2,0.2, 200)
norm = colors.BoundaryNorm(bounds, cmap.N)
flag = True
def data(t):
    # Create the mesh in polar coordinates and compute corresponding Z.
    Z = res1[:,t]*R
    ax.clear()
    # Express the mesh in the cartesian system.
    
    surf = ax.plot_surface(X, Y, Z,cmap=cmap, norm=norm ,lw=3)
    if(t==1):
        fig.colorbar(surf, shrink=0.5, aspect=10)
    # ax.set_xlim(0, 4)
    ax.set_zlim(0, 1)
    # ax.set_ylim(0, 2)
    

fig = plt.figure()
ax = fig.gca(projection='3d')

r=np.linspace(0, Test.x, num = Test.len_X)
p = np.linspace(0,2*np.pi, num = Test.len_X)
R, P = np.meshgrid(r, p)
X, Y = R*np.cos(P), R*np.sin(P)

ani = FuncAnimation(fig, data, frames=Frames, interval=20  )
 
 
ani.save('mesh.gif', writer='imagemagick')

###########################






