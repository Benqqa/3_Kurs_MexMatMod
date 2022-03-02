# -*- coding: utf-8 -*-
"""
Created on Tue Nov  2 10:26:37 2021

@author: maxbo
"""
import numpy as np
import sympy as sym
import matplotlib.pyplot as plt
from matplotlib import cm
from mpl_toolkits.mplot3d import Axes3D
class Tochka():
    def __init__(self,x,y,z):
        self.x=x
        self.y=y
        self.z=z
prisma={}
prisma[0]=Tochka(0,0,0)
prisma[1]=Tochka(0,1,0)
prisma[2]=Tochka(1,0,0)
prisma[3]=Tochka(0,0,1)
prisma[4]=Tochka(1,0,1)
prisma[5]=Tochka(0,1,1)
print(prisma[5].z)
X=np.zeros((len(prisma),len(prisma)))
X1 = sym.Matrix( len(prisma), len(prisma),range(len(prisma)**2))

k=0
for t in prisma:
    X1[k,0]=1
    X1[k,1]=prisma[t].x
    X1[k,2]=prisma[t].y
    X1[k,3]=prisma[t].z
    X1[k,4]=prisma[t].x*prisma[t].z
    X1[k,5]=prisma[t].y*prisma[t].z
    k+=1
print(X1)
A1=X1.inv()
print(A1)
x,y,z=sym.symbols('x y z')
p1=sym.Matrix([1,x,y,z,x*z,y*z])
print(p1)
print("aaaaa",p1.T)
P1 = sym.Matrix( len(prisma), len(prisma),range(len(prisma)**2))
for i in prisma:
    P1[i,:]=[[1,x,y,z,x*z,y*z]]
print(P1)
res=A1*p1
print("rez",res)
fig, ax = plt.subplots()
ax.plot(res)
plt.show()
# X=np.linspace(0, 1, num = Test.len_X)
# Y=np.linspace(0, 1, num = Test.len_T)
# X,Y = np.meshgrid(X,Y)
# fig, ax = plt.subplots(subplot_kw={"projection": "3d"})
# surf = ax.plot_surface(X, Y, res1, cmap='coolwarm',linewidth=0)
# X[:,0]=1
# k=0
# for t in prisma:
#     X[k,0]=1
#     X[k,1]=prisma[t].x
#     X[k,2]=prisma[t].y
#     X[k,3]=prisma[t].z
#     X[k,4]=prisma[t].x*prisma[t].z
#     X[k,5]=prisma[t].y*prisma[t].z
#     k+=1
# print(X)
# A=np.linalg.inv(X)
# print(A)
# x=sym.symbols('x')
# y=sym.symbols('y')
# z=sym.symbols('z')
# P=sym.Array([1,x,y,z,x*z,y*z])
# A1=sym.Array(range(len(prisma)**2), (len(prisma), len(prisma)))
# print(A1.tomatrix())
# # P[0]=1
# # P[1]=x
# # P[2]=y
# # P[3]=z
# # P[4]=x*z
# # P[5]=y*z
# for i in range(len(prisma)):
#     for j in range(len(prisma)):
#         A1[i,j]=X[i,j]

# print(A1)